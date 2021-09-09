/*
    Copyright (C) 2016 Tomas Flouri, Kassian Kobert

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

static void fill_parent_scaler(unsigned int scaler_size,
                               unsigned int * parent_scaler,
                               const unsigned int * left_scaler,
                               const unsigned int * right_scaler)
{
  unsigned int i;

  if (!left_scaler && !right_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * scaler_size);
  else if (left_scaler && right_scaler)
  {
    memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    for (i = 0; i < scaler_size; ++i)
      parent_scaler[i] += right_scaler[i];
  }
  else
  {
    if (left_scaler)
      memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    else
      memcpy(parent_scaler, right_scaler, sizeof(unsigned int) * scaler_size);
  }
}


PLL_EXPORT void pll_core_update_partial_ti_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                double * parent_clv,
                                                unsigned int * parent_scaler,
                                                const unsigned char * left_tipchars,
                                                const double * right_clv,
                                                const double * left_matrix,
                                                const double * right_matrix,
                                                const unsigned int * right_scaler,
                                                const pll_state_t * tipmap,
                                                unsigned int tipmap_size,
                                                unsigned int attrib)
{
  unsigned int i,j,k,n;

  const double * lmat;
  const double * rmat;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = states_padded * rate_cats;

  pll_state_t lstate;

  /* dedicated functions for 4x4 matrices (DNA) */
  if (states == 4)
  {
    /* no AVX2 kernel so far; rollback to AVX */
    pll_core_update_partial_ti_4x4_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_tipchars,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       right_scaler,
                                       attrib);
    return;
  }

  /* dedicated functions for 20x20 matrices (AA) */
  if (states == 20)
  {
    pll_core_update_partial_ti_20x20_avx2(sites,
                                          rate_cats,
                                          parent_clv,
                                          parent_scaler,
                                          left_tipchars,
                                          right_clv,
                                          left_matrix,
                                          right_matrix,
                                          right_scaler,
                                          tipmap,
                                          tipmap_size,
                                          attrib);
    return;
  }

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, NULL, right_scaler);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  __m256i mask;

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;

    scale_mask = init_mask;

    lstate = tipmap[left_tipchars[n]];

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {

        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;
        __m256d v_rclv;

        /* point to the four rows of the left matrix */
        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* set position of least significant bit in character state */
        register int lsb = 0;

        /* iterate over quadruples of columns */
        for (j = 0; j < states_padded; j += 4)
        {

          /* set mask */
          mask = _mm256_set_epi64x(
                    ((lstate >> (lsb+3)) & 1) ? ~0 : 0,
                    ((lstate >> (lsb+2)) & 1) ? ~0 : 0,
                    ((lstate >> (lsb+1)) & 1) ? ~0 : 0,
                    ((lstate >> (lsb+0)) & 1) ? ~0 : 0);

          if ((lstate >> lsb) & 0b1111)
          {
            v_mat    = _mm256_maskload_pd(lm0,mask);
            v_terma0 = _mm256_add_pd(v_terma0,v_mat);

            v_mat    = _mm256_maskload_pd(lm1,mask);
            v_terma1 = _mm256_add_pd(v_terma1,v_mat);

            v_mat    = _mm256_maskload_pd(lm2,mask);
            v_terma2 = _mm256_add_pd(v_terma2,v_mat);

            v_mat    = _mm256_maskload_pd(lm3,mask);
            v_terma3 = _mm256_add_pd(v_terma3,v_mat);
          }

          lsb += 4;

          lm0 += 4;
          lm1 += 4;
          lm2 += 4;
          lm3 += 4;

          v_rclv    = _mm256_load_pd(right_clv+j);

          /* row 0 */
          v_mat    = _mm256_load_pd(rm0);
          v_termb0 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb0);
          rm0 += 4;

          /* row 1 */
          v_mat    = _mm256_load_pd(rm1);
          v_termb1 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb1);
          rm1 += 4;

          /* row 2 */
          v_mat    = _mm256_load_pd(rm2);
          v_termb2 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb2);
          rm2 += 4;

          /* row 3 */
          v_mat    = _mm256_load_pd(rm3);
          v_termb3 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb3);
          rm3 += 4;
        }

        /* point pmatrix to the next four rows */ 
        lmat = lm3;
        rmat = rm3;

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);

      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      parent_clv += states_padded;
      right_clv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

#define COMPUTE_TI_QCOL(q, offset) \
/* row 0 */ \
v_mat    = _mm256_load_pd(rm0 + offset); \
v_termb0 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb0); \
 \
/* row 1 */ \
v_mat    = _mm256_load_pd(rm1 + offset); \
v_termb1 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb1); \
\
/* row 2 */ \
v_mat    = _mm256_load_pd(rm2 + offset); \
v_termb2 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb2); \
\
/* row 3 */ \
v_mat    = _mm256_load_pd(rm3 + offset); \
v_termb3 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb3);

PLL_EXPORT
void pll_core_update_partial_ti_20x20_avx2(unsigned int sites,
                                           unsigned int rate_cats,
                                           double * parent_clv,
                                           unsigned int * parent_scaler,
                                           const unsigned char * left_tipchar,
                                           const double * right_clv,
                                           const double * left_matrix,
                                           const double * right_matrix,
                                           const unsigned int * right_scaler,
                                           const pll_state_t * tipmap,
                                           unsigned int tipmap_size,
                                           unsigned int attrib)
{
  unsigned int states = 20;
  unsigned int states_padded = states;
  unsigned int maxstates = tipmap_size;
  unsigned int i,j,k,n,m;

  const double * lmat;
  const double * rmat;

  unsigned int span_padded = states_padded * rate_cats;
  unsigned int lstate;

  __m256d xmm0,xmm1,xmm2,xmm3;

  /* precompute a lookup table of four values per entry (one for each state),
     for all 16 states (including ambiguities) and for each rate category. */
  double * lookup = pll_aligned_alloc(maxstates*span_padded*sizeof(double),
                                      PLL_ALIGNMENT_AVX);
  if (!lookup)
  {
    /* TODO: in the highly unlikely event that allocation fails, we should
       resort to a non-lookup-precomputation version of this function,
       available at commit e.g.  a4fc873fdc65741e402cdc1c59919375143d97d1 */
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate space for precomputation.");
    return;
  }

  double * ptr = lookup;

  /* precompute left-side values and store them in lookup table */
  for (j = 0; j < maxstates; ++j)
  {
    lmat = left_matrix;

    // we can do this since the number of state is fixed to 20 (and thus < 32)
    unsigned int state = (unsigned int) tipmap[j];

    unsigned int scount = PLL_POPCNT32(state);

    if (scount == 1)
    {
      /* special case for non-ambiguous states */
      int ss = PLL_CTZ32(state);

      for (n = 0; n < rate_cats; ++n)
      {
        for (i = 0; i < states; ++i)
        {
          ptr[i] = lmat[ss];

          lmat += states;
        }

        ptr += states;
      }
    }
    else if (scount == states)
    {
      /* special case for gaps */
      for (n = 0; n < rate_cats; ++n)
      {
        for (i = 0; i < states; ++i)
        {
          double terml = 0;
          for (m = 0; m < states; ++m)
          {
            terml += lmat[m];
          }
          ptr[i] = terml;

          lmat += states;
        }

        ptr += states;
      }
    }
    else
    {
      /* generic case for ambiguous states */
      for (n = 0; n < rate_cats; ++n)
      {
        for (i = 0; i < states; ++i)
        {
          double terml = 0;
          for (m = 0; m < states; ++m)
          {
            if ((state>>m) & 1)
            {
              terml += lmat[m];
            }
          }
          ptr[i] = terml;

          lmat += states;
        }

        ptr += states;
      }
    }
  }

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, NULL, right_scaler);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  /* iterate over sites and compute CLV entries */
  for (n = 0; n < sites; ++n)
  {
    rmat = right_matrix;

    scale_mask = init_mask;

    lstate = (unsigned int) left_tipchar[n];

    double * state_lookup = lookup + lstate*span_padded;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      __m256d v_rclv[5];
      for (i = 0; i < 5; ++i)
      {
        v_rclv[i]    = _mm256_load_pd(right_clv + i*4);
      }

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* iterate over quadruples of columns */
//        COMPUTE_TI_QCOL_FIRST;
        COMPUTE_TI_QCOL(0, 0);
        COMPUTE_TI_QCOL(1, 4);
        COMPUTE_TI_QCOL(2, 8);
        COMPUTE_TI_QCOL(3, 12);
        COMPUTE_TI_QCOL(4, 16);

        /* point pmatrix to the next four rows */
        rmat = rm3 + states_padded;

        /* load x from precomputed lookup table */
        __m256d v_terma_sum = _mm256_load_pd(state_lookup+i);

        /* compute termb */
        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);
      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      rmat -= displacement;

      parent_clv += states_padded;
      right_clv  += states_padded;
      state_lookup += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
  pll_aligned_free(lookup);
}

#define COMPUTE_II_QCOL(q, offset) \
/* row 0 */ \
v_mat    = _mm256_load_pd(lm0 + offset); \
v_terma0 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_terma0); \
v_mat    = _mm256_load_pd(rm0 + offset); \
v_termb0 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb0); \
 \
/* row 1 */ \
v_mat    = _mm256_load_pd(lm1 + offset); \
v_terma1 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_terma1); \
v_mat    = _mm256_load_pd(rm1 + offset); \
v_termb1 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb1); \
\
/* row 2 */ \
v_mat    = _mm256_load_pd(lm2 + offset); \
v_terma2 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_terma2); \
v_mat    = _mm256_load_pd(rm2 + offset); \
v_termb2 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb2); \
\
/* row 3 */ \
v_mat    = _mm256_load_pd(lm3 + offset); \
v_terma3 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_terma3); \
v_mat    = _mm256_load_pd(rm3 + offset); \
v_termb3 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_termb3);


PLL_EXPORT void pll_core_update_partial_ii_20x20_avx2(unsigned int sites,
                                                unsigned int rate_cats,
                                                double * parent_clv,
                                                unsigned int * parent_scaler,
                                                const double * left_clv,
                                                const double * right_clv,
                                                const double * left_matrix,
                                                const double * right_matrix,
                                                const unsigned int * left_scaler,
                                                const unsigned int * right_scaler,
                                                unsigned int attrib)
{
  unsigned int i,k,n;

  const double * lmat;
  const double * rmat;

  unsigned int states = 20;
  unsigned int states_padded = states;
  unsigned int span_padded = states_padded * rate_cats;

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, left_scaler, right_scaler);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      __m256d v_lclv[5];
      __m256d v_rclv[5];
      for (i = 0; i < 5; ++i)
      {
        v_lclv[i]    = _mm256_load_pd(left_clv + i*4);
        v_rclv[i]    = _mm256_load_pd(right_clv + i*4);
      }

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_mat;

        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        /* point to the four rows of the left matrix */
        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* iterate over quadruples of columns */
        COMPUTE_II_QCOL(0, 0);
        COMPUTE_II_QCOL(1, 4);
        COMPUTE_II_QCOL(2, 8);
        COMPUTE_II_QCOL(3, 12);
        COMPUTE_II_QCOL(4, 16);

        /* point pmatrix to the next four rows */
        lmat = lm3 + states_padded;
        rmat = rm3 + states_padded;

        /* compute terma */

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);
      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      parent_clv += states_padded;
      left_clv   += states_padded;
      right_clv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_core_update_partial_ii_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                double * parent_clv,
                                                unsigned int * parent_scaler,
                                                const double * left_clv,
                                                const double * right_clv,
                                                const double * left_matrix,
                                                const double * right_matrix,
                                                const unsigned int * left_scaler,
                                                const unsigned int * right_scaler,
                                                unsigned int attrib)
{
  unsigned int i,j,k,n;

  const double * lmat;
  const double * rmat;

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = states_padded * rate_cats;

  /* dedicated functions for 4x4 and 20x20 matrices */
  if (states == 4)
  {
    /* TODO: Implement avx2 4x4 case */
    pll_core_update_partial_ii_4x4_avx(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_clv,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       left_scaler,
                                       right_scaler,
                                       attrib);
    return;
  }
  else if (states == 20)
  {
    pll_core_update_partial_ii_20x20_avx2(sites,
                                       rate_cats,
                                       parent_clv,
                                       parent_scaler,
                                       left_clv,
                                       right_clv,
                                       left_matrix,
                                       right_matrix,
                                       left_scaler,
                                       right_scaler,
                                       attrib);
    return;
  }

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    const size_t scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;
    /* add up the scale vector of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, left_scaler, right_scaler);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;
        __m256d v_lclv;
        __m256d v_rclv;

        /* point to the four rows of the left matrix */
        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* iterate over quadruples of columns */
        for (j = 0; j < states_padded; j += 4)
        {
          v_lclv    = _mm256_load_pd(left_clv+j);
          v_rclv    = _mm256_load_pd(right_clv+j);

          /* row 0 */
          v_mat    = _mm256_load_pd(lm0);
          v_terma0 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma0);

          v_mat    = _mm256_load_pd(rm0);
          v_termb0 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb0);
          lm0 += 4;
          rm0 += 4;

          /* row 1 */
          v_mat    = _mm256_load_pd(lm1);
          v_terma1 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma1);

          v_mat    = _mm256_load_pd(rm1);
          v_termb1 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb1);
          lm1 += 4;
          rm1 += 4;

          /* row 2 */
          v_mat    = _mm256_load_pd(lm2);
          v_terma2 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma2);

          v_mat    = _mm256_load_pd(rm2);
          v_termb2 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb2);
          lm2 += 4;
          rm2 += 4;

          /* row 3 */
          v_mat    = _mm256_load_pd(lm3);
          v_terma3 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma3);

          v_mat    = _mm256_load_pd(rm3);
          v_termb3 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb3);

          lm3 += 4;
          rm3 += 4;
        }

        /* point pmatrix to the next four rows */ 
        lmat = lm3;
        rmat = rm3;

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);
      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      parent_clv += states_padded;
      left_clv   += states_padded;
      right_clv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

PLL_EXPORT void pll_core_update_partial_repeats_20x20_avx2(unsigned int parent_sites,
                                                           unsigned int left_sites,
                                                           unsigned int right_sites,
                                                           unsigned int rate_cats,
                                                           double * parent_clv,
                                                           unsigned int * parent_scaler,
                                                           const double * left_clv,
                                                           const double * right_clv,
                                                           const double * left_matrix,
                                                           const double * right_matrix,
                                                           const unsigned int * left_scaler,
                                                           const unsigned int * right_scaler,
                                                           const unsigned int * parent_id_site,
                                                           const unsigned int * left_site_id,
                                                           const unsigned int * right_site_id,
                                                           double * bclv_buffer,
                                                           unsigned int attrib)
{
  unsigned int i,k,n;

  const double * lmat;
  const double * rmat;

  unsigned int states = 20;
  unsigned int states_padded = states;
  unsigned int span_padded = states_padded * rate_cats;

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    /* add up the scale vector of the two children if available */
    if (scale_mode == 2)
      pll_fill_parent_scaler_repeats_per_rate(parent_sites, rate_cats, parent_scaler, parent_id_site,
        left_scaler, left_site_id, right_scaler, right_site_id);
    else
      pll_fill_parent_scaler_repeats(parent_sites, parent_scaler, parent_id_site,
        left_scaler, left_site_id, right_scaler, right_site_id);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  /* compute CLV */
  for (n = 0; n < parent_sites; ++n)
  {
    unsigned int site = PLL_GET_SITE(parent_id_site, n);
    unsigned int lid = PLL_GET_ID(left_site_id, site);
    unsigned int rid = PLL_GET_ID(right_site_id, site);
    const double *lclv = &left_clv[lid * span_padded];
    const double *rclv = &right_clv[rid * span_padded];
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      __m256d v_lclv[5];
      __m256d v_rclv[5];
      for (i = 0; i < 5; ++i)
      {
        v_lclv[i]    = _mm256_load_pd(lclv + i*4);
        v_rclv[i]    = _mm256_load_pd(rclv + i*4);
      }

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;

        /* point to the four rows of the left matrix */
        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* iterate over quadruples of columns */
        COMPUTE_II_QCOL(0, 0);
        COMPUTE_II_QCOL(1, 4);
        COMPUTE_II_QCOL(2, 8);
        COMPUTE_II_QCOL(3, 12);
        COMPUTE_II_QCOL(4, 16);

        /* point pmatrix to the next four rows */
        lmat = lm3 + states_padded;
        rmat = rm3 + states_padded;

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);
      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      parent_clv += states_padded;
      lclv   += states_padded;
      rclv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}


PLL_EXPORT void pll_core_update_partial_repeats_generic_avx2(unsigned int states,
                                                             unsigned int parent_sites,
                                                             unsigned int left_sites,
                                                             unsigned int right_sites,
                                                             unsigned int rate_cats,
                                                             double * parent_clv,
                                                             unsigned int * parent_scaler,
                                                             const double * left_clv,
                                                             const double * right_clv,
                                                             const double * left_matrix,
                                                             const double * right_matrix,
                                                             const unsigned int * left_scaler,
                                                             const unsigned int * right_scaler,
                                                             const unsigned int * parent_id_site,
                                                             const unsigned int * left_site_id,
                                                             const unsigned int * right_site_id,
                                                             double * bclv_buffer,
                                                             unsigned int attrib)
{
  unsigned int i,j,k,n;

  const double * lmat;
  const double * rmat;

  if (states == 20)
  {
    return pll_core_update_partial_repeats_20x20_avx2(parent_sites,
                                                      left_sites,
                                                      right_sites,
                                                      rate_cats,
                                                      parent_clv,
                                                      parent_scaler,
                                                      left_clv,
                                                      right_clv,
                                                      left_matrix,
                                                      right_matrix,
                                                      left_scaler,
                                                      right_scaler,
                                                      parent_id_site,
                                                      left_site_id,
                                                      right_site_id,
                                                      bclv_buffer,
                                                      attrib);
  }

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = states_padded * rate_cats;

  /* scaling-related stuff */
  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int scale_mask;
  unsigned int init_mask;
  __m256d v_scale_threshold = _mm256_set1_pd(PLL_SCALE_THRESHOLD);
  __m256d v_scale_factor = _mm256_set1_pd(PLL_SCALE_FACTOR);

  if (!parent_scaler)
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  else
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 0xF : 0;
    /* add up the scale vector of the two children if available */
    if (scale_mode == 2) 
      pll_fill_parent_scaler_repeats_per_rate(parent_sites, rate_cats, parent_scaler, parent_id_site, 
        left_scaler, left_site_id, right_scaler, right_site_id);
    else
      pll_fill_parent_scaler_repeats(parent_sites, parent_scaler, parent_id_site, 
        left_scaler, left_site_id, right_scaler, right_site_id);
  }

  size_t displacement = (states_padded - states) * (states_padded);

  /* compute CLV */
  for (n = 0; n < parent_sites; ++n)
  {
    unsigned int site = PLL_GET_SITE(parent_id_site, n);
    unsigned int lid = PLL_GET_ID(left_site_id, site);
    unsigned int rid = PLL_GET_ID(right_site_id, site);
    const double *lclv = &left_clv[lid * span_padded];
    const double *rclv = &right_clv[rid * span_padded];
    lmat = left_matrix;
    rmat = right_matrix;
    scale_mask = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_mask = 0xF;

      /* iterate over quadruples of rows */
      for (i = 0; i < states_padded; i += 4)
      {
        __m256d v_terma0 = _mm256_setzero_pd();
        __m256d v_termb0 = _mm256_setzero_pd();
        __m256d v_terma1 = _mm256_setzero_pd();
        __m256d v_termb1 = _mm256_setzero_pd();
        __m256d v_terma2 = _mm256_setzero_pd();
        __m256d v_termb2 = _mm256_setzero_pd();
        __m256d v_terma3 = _mm256_setzero_pd();
        __m256d v_termb3 = _mm256_setzero_pd();

        __m256d v_mat;
        __m256d v_lclv;
        __m256d v_rclv;

        /* point to the four rows of the left matrix */
        const double * lm0 = lmat;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;

        /* point to the four rows of the right matrix */
        const double * rm0 = rmat;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;

        /* iterate over quadruples of columns */
        for (j = 0; j < states_padded; j += 4)
        {
          v_lclv    = _mm256_load_pd(lclv+j);
          v_rclv    = _mm256_load_pd(rclv+j);

          /* row 0 */
          v_mat    = _mm256_load_pd(lm0);
          v_terma0 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma0);

          v_mat    = _mm256_load_pd(rm0);
          v_termb0 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb0);
          lm0 += 4;
          rm0 += 4;

          /* row 1 */
          v_mat    = _mm256_load_pd(lm1);
          v_terma1 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma1);

          v_mat    = _mm256_load_pd(rm1);
          v_termb1 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb1);
          lm1 += 4;
          rm1 += 4;

          /* row 2 */
          v_mat    = _mm256_load_pd(lm2);
          v_terma2 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma2);

          v_mat    = _mm256_load_pd(rm2);
          v_termb2 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb2);
          lm2 += 4;
          rm2 += 4;

          /* row 3 */
          v_mat    = _mm256_load_pd(lm3);
          v_terma3 = _mm256_fmadd_pd(v_mat, v_lclv, v_terma3);

          v_mat    = _mm256_load_pd(rm3);
          v_termb3 = _mm256_fmadd_pd(v_mat, v_rclv, v_termb3);

          lm3 += 4;
          rm3 += 4;
        }

        /* point pmatrix to the next four rows */ 
        lmat = lm3;
        rmat = rm3;

        __m256d xmm0 = _mm256_unpackhi_pd(v_terma0,v_terma1);
        __m256d xmm1 = _mm256_unpacklo_pd(v_terma0,v_terma1);

        __m256d xmm2 = _mm256_unpackhi_pd(v_terma2,v_terma3);
        __m256d xmm3 = _mm256_unpacklo_pd(v_terma2,v_terma3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_terma_sum = _mm256_add_pd(xmm2,xmm3);

        /* compute termb */

        xmm0 = _mm256_unpackhi_pd(v_termb0,v_termb1);
        xmm1 = _mm256_unpacklo_pd(v_termb0,v_termb1);

        xmm2 = _mm256_unpackhi_pd(v_termb2,v_termb3);
        xmm3 = _mm256_unpacklo_pd(v_termb2,v_termb3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_termb_sum = _mm256_add_pd(xmm2,xmm3);

        __m256d v_prod = _mm256_mul_pd(v_terma_sum,v_termb_sum);

        /* check if scaling is needed for the current rate category */
        __m256d v_cmp = _mm256_cmp_pd(v_prod, v_scale_threshold, _CMP_LT_OS);
        rate_mask = rate_mask & _mm256_movemask_pd(v_cmp);

        _mm256_store_pd(parent_clv+i, v_prod);
      }

      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_mask == 0xF)
        {
          for (i = 0; i < states_padded; i += 4)
          {
            __m256d v_prod = _mm256_load_pd(parent_clv + i);
            v_prod = _mm256_mul_pd(v_prod, v_scale_factor);
            _mm256_store_pd(parent_clv + i, v_prod);
          }
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        scale_mask = scale_mask & rate_mask;

      /* reset pointers to point to the start of the next p-matrix, as the
         vectorization assumes a square states_padded * states_padded matrix,
         even though the real matrix is states * states_padded */
      lmat -= displacement;
      rmat -= displacement;

      parent_clv += states_padded;
      lclv   += states_padded;
      rclv  += states_padded;
    }

    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scale_mask == 0xF)
    {
      parent_clv -= span_padded;
      for (i = 0; i < span_padded; i += 4)
      {
        __m256d v_prod = _mm256_load_pd(parent_clv + i);
        v_prod = _mm256_mul_pd(v_prod,v_scale_factor);
        _mm256_store_pd(parent_clv + i, v_prod);
      }
      parent_clv += span_padded;
      parent_scaler[n] += 1;
    }
  }
}

