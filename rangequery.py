from __future__ import print_function

import dataanalysis as da
import ddosa

import os
import re

import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy import units as u

import pilton

def converttime(informat,intime,outformat):
    ct=pilton.heatool("converttime")
    ct['informat']=informat
    ct['intime']=intime
    ct['outformat']=outformat
    ct.run()

    c=dict(re.findall("Output Time\((.*)\): (.*)",ct.output))
    if outformat=="":
        return c
    else:
        return c[outformat]

class TimeDirectionScWList(ddosa.DataAnalysis):
    coordinates=dict(RA=None,DEC=None,radius=None)
    timespan=dict(T1=None,T2=None)
    max_pointings=10

    allow_alias=True

    def get_version(self):
        v=self.get_signature()+"."+self.version

        try:
            v+="c_%(RA).5lg_%(DEC).5lg_%(radius).5lg.t_%(T1)s_%(T2)s_maxpoint%(max_pointings)i"%(
                            dict(self.coordinates.items()+ \
                                 self.timespan.items()+ \
                                 dict(max_pointings=(self.max_pointings if self.max_pointings is not None else 1000)).items())
                        )
        except Exception:
            print("failed coordinates",self.coordinates)
            print("failed timespan", self.timespan)
            print("failed max_pointings", self.max_pointings)
            raise

        return v

    def main(self):
        scw_index=fits.open(os.environ['INTEGRAL_DATA']+"/idx/scw/GNRL-SCWG-GRP-IDX.fits")[1].data

        scx=SkyCoord(scw_index['RA_SCX'],scw_index['DEC_SCX'],unit="deg")

        target=SkyCoord(self.coordinates['RA'],self.coordinates['DEC'],unit="deg")

        m_avail=scw_index['SW_TYPE']=="POINTING"

        m_coord=scx.separation(target).deg<self.coordinates['radius']

        print("for coordinates",target,"found",sum(m_coord & m_avail),"within",self.coordinates['radius'],"degrees")
        
        t1_ijd=float(converttime("UTC",self.timespan['T1'],"IJD"))
        t2_ijd=float(converttime("UTC",self.timespan['T2'],"IJD"))

        t1=scw_index['TSTART']
        t2=scw_index['TSTOP']
        m_time=(t1>t1_ijd) & (t2<t2_ijd)

        print("for time",t1_ijd,t2_ijd,"found",sum(m_coord & m_avail & m_time), "total range was",t1.min(),t2.max(),"(",(t2.max()-t1.min())/365,"years)")
            
        pre_selection=scw_index['SWID'][m_coord & m_avail & m_time]
        selection=[]

        for scwid in pre_selection:
            evtsfn=os.environ['INTEGRAL_DATA']+"/scw/%s/%s.001/isgri_events.fits.gz"%(scwid[:4],scwid)
            print("searching for",evtsfn)
            if not os.path.exists(evtsfn):
                print("skipping",scwid)
                continue

            selection.append(scwid)

            if self.max_pointings is not None and len(selection)>=self.max_pointings:
                print("choosing only first",self.max_pointings)
                break

        print("selection:",selection)

        self.scwlistdata=[ddosa.ScWData(input_scwid=scwid+".001") for scwid in selection]
