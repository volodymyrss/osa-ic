from __future__ import print_function

import astropy.io.fits as pyfits
import argparse
import pilton
from ddosa import remove_withtemplate
import os
import numpy as np


class ICTree(object):
    def __init__(self):
        self.icroot=os.environ['CURRENT_IC']

    @property
    def ibisicroot(self):
        return self.icroot+"/ic/ibis/mod"

    @property
    def idxicroot(self):
        return self.icroot+"/idx/ic/"

    @property
    def icmaster(self):
        return self.get_icmaster()
    
    def get_icmaster(self,version=None):
        return self.icroot+"/idx/ic/ic_master_file"+("_"+version if version is not None else "")+".fits"

    def DS_to_fn_prefix(self,DS):
        return DS.lower().replace("-","_")

    def DS_to_fn(self,DS,serial=0):
        return self.ibisicroot+"/"+self.DS_to_fn_prefix(DS)+"_%.4i.fits"%serial

    def DS_to_idx_fn(self,DS):
        return self.idxicroot+"/"+DS+"-IDX.fits"

    def create_index_empty(self,DS):
        dc=pilton.heatool("dal_create")
        dc["obj_name"]=self.DS_to_idx_fn(DS)
        dc["template"]=DS+"-IDX.tpl"
        remove_withtemplate(dc["obj_name"].value+"("+dc["template"].value+")")
        dc.run()
        
    def get_file_DS(self,fn):
        f=pyfits.open(fn)
        if len(f)>2:
            print("we refuse to deal with indexed IC ds, as they increase the amount of suffering in the world")
            raise Exception("too many extensions %i %s"%(len(f),repr(f)))
        if len(f)<2:
            print("")
            raise Exception("too few extensions %i %s"%(len(f),repr(f)))
        return f[1].header['EXTNAME']


    def attach_ds(self,fn,serial=0):
        DS=self.get_file_DS(fn)
        f_ds=pyfits.open(fn)
        f_ds.writeto(self.DS_to_fn(DS,serial),clobber=True)

        da=pilton.heatool("dal_attach")
        da['Parent']=self.DS_to_idx_fn(DS)
        da['Child1']=self.DS_to_fn(DS,serial)
        da.run()

        f_idx=pyfits.open(da['Parent'].value)

        for k in ['VERSION','VSTART','VSTOP']:
            f_idx[1].data[-1][k]=f_ds[1].header[k]
            print(k,f_idx[1].data[-1][k],f_ds[1].header[k])

        f_idx.writeto(da['Parent'].value,clobber=True)

        dv=pilton.heatool("dal_verify")
        dv["indol"]=self.DS_to_idx_fn(DS)
        dv['checksums']="yes"
        dv['backpointers']="yes"
        dv['detachother']="yes"
        dv.run()

    def init_icmaster(self):
        f=pyfits.open(self.get_icmaster("osa102"))
        
        print("columns was:",len(f[3].columns))

        nhdu=pyfits.BinTableHDU.from_columns(f[3].columns+pyfits.ColDefs([pyfits.Column("ISGR_L2RE_MOD","1I"),pyfits.Column("ISGR_MCEC_MOD","1I"),pyfits.Column("ISGR_EFFC_MOD","1I")]))
        #nhdu=pyfits.BinTableHDU.from_columns(f[3].columns+pyfits.ColDefs([pyfits.Column("ISGR_LUT2_MOD","1I"),pyfits.Column("ISGR_L2RE_MOD","1I"),pyfits.Column("ISGR_MCEC_MOD","1I"),pyfits.Column("ISGR_EFFC_MOD","1I")]))

        for k,v in f[3].header.items():
            if k.startswith("TFORM"): continue
            if k.startswith("TTYPE"): continue
            if k.startswith("TFIELDS"): continue
            if k.startswith("NAXIS"): continue
            #print ":",k,"x",v
            nhdu.header[k]=v

        nhdu.data[0]['ISGR_RISE_MOD']=2
        nhdu.data[0]['ISGR_L2RE_MOD']=1
        nhdu.data[0]['ISGR_MCEC_MOD']=1
        nhdu.data[0]['ISGR_EFFC_MOD']=1

        f[3]=nhdu

        print("columns now:",len(f[3].columns))

        f.writeto(self.get_icmaster(),clobber=True)

    def create_index_from_list(self,fns=None,fns_list=None):
        if fns_list is None:
            if fns is None:
                raise Exception("what?")

            fns_list=tempfile.mkstemp()
            open(fns_list,"w").write("\n".join(fns))

        else:
            if fns is not None:
                raise Exception("what?")


        for rev in [239,665,1516]:
            infn="/sps/integral/data/reduced/ddcache//byrev/%.4i/ISGRI_RISE_MOD.v0//2e84304c/isgr_rise_mod_%.4i.fits.gz"%(rev,rev)
            f=pyfits.open(infn)

            outfn="/sps/integral/data/ic/ic_snapshot_20140321/ic/ibis/mod/"+f[1].header['FILENAME']

            f.writeto(outfn,clobber=True)

            listf.write(outfn+"\n")


        listf.close()

        da=pilton.heatool("txt2idx")
        da['index']="/sps/integral/data/ic/ic_snapshot_20140321/idx/ic/ISGR-RISE-MOD-IDX.fits"
        da['template']="ISGR-RISE-MOD-IDX.tpl"
        da['update']=1
        da['element']=listfn
        da.run()

        f=pyfits.open(da['index'].value)
        f[1].header['CREATOR']=""
        f[1].header['CONFIGUR']=""
        f.writeto(da['index'].value,clobber=True)

        dv=pilton.heatool("dal_verify")

    def attach_to_master(self,DS):
        da=pilton.heatool("dal_attach")
        da['Parent']=self.icmaster
        da['Child1']=self.DS_to_idx_fn(DS)
        da.run()


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('icmembers', metavar='new IC file', nargs='+',
            help='new IC file')
    #parser.add_argument('clobber-index', action='store_true',
    #        help='clobber index')

    args = parser.parse_args()

    ictree=ICTree()
    ictree.init_icmaster()

    for icmember in args.icmembers:
        DS=ictree.get_file_DS(icmember)
        print(icmember,"as",DS)

        ictree.create_index_empty(DS)
        ictree.attach_ds(icmember)
        ictree.attach_to_master(DS)

if __name__ == "__main__":
    main()


