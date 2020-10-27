#!/bin/bash

rev=${1:-*}
shift 1

osa-ic \
    -f <(ls osa11_deployment/deployment/ic/ic/ibis/rsp/isgr_arf_rsp_0052.fits) \
    -f <(ls -tr ddcache/byrev/${rev}/ISGR_RISE_MOD_Revolution.*/*/isgr_rise_mod_*.fits.gz | tail -1) \
    -f <(ls -tr ddcache/byrev/${rev}/ISGR_EFFC_MOD_Revolution.*/ResponseRev.v3_r1_1_r2_0/*/isgr_effc_mod_*.fits.gz | tail -1) \
    -f <(ls -tr ddcache/byrev/${rev}/ISGR_MCEC_MOD_Revolution.v0/33f56695/isgr_mcec_mod_*.fits.gz | tail -1) \
    -f <(ls -tr ddcache/byrev/${rev}/ISGR_L2RE_MOD_Revolution.v0/33f56695/isgr_l2re_mod_*.fits.gz | tail -1) \
    -f <(ls -tr ddcache/byrev/${rev}/ResponseIC_Revolution.v0/ResponseRev.v3_r1_1_r2_0/c7d24cae/isgr_rmf_rsp_*.fits.gz | tail -1) \
    $@
    
