"""Reduced hydraulic and electrical-heating scenarios, not CFD or thermal validation."""
from pathlib import Path
import argparse,csv,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]

def hydraulic(flow_uL_min=450.,length_m=.06,height_m=.0003,width_m=.012,
              density_kg_m3=786.,viscosity_Pa_s=.00035,diffusivity_m2_s=1.1e-9):
    values=[flow_uL_min,length_m,height_m,width_m,density_kg_m3,viscosity_Pa_s,diffusivity_m2_s]
    if any(not math.isfinite(x) or x<=0 for x in values):raise ValueError('Positive finite physical inputs required')
    q=flow_uL_min*1e-9/60;area=height_m*width_m;u=q/area;volume=area*length_m;diameter=2*area/(height_m+width_m)
    dp=12*viscosity_Pa_s*length_m*q/(width_m*height_m**3)
    return dict(flow_uL_min=flow_uL_min,length_m=length_m,height_m=height_m,width_m=width_m,density_kg_m3=density_kg_m3,viscosity_Pa_s=viscosity_Pa_s,diffusivity_m2_s=diffusivity_m2_s,flow_m3_s=q,volume_mL=volume*1e6,velocity_m_s=u,hydraulic_diameter_m=diameter,residence_time_s=volume/q,Re_hydraulic=density_kg_m3*u*diameter/viscosity_Pa_s,Pe_height=u*height_m/diffusivity_m2_s,Pe_length=u*length_m/diffusivity_m2_s,transverse_diffusion_scale_s=height_m**2/diffusivity_m2_s,parallel_plate_pressure_drop_Pa=dp,ideal_hydraulic_power_W=dp*q,aspect_ratio_H_over_W=height_m/width_m,evidence_role='uncalibrated_hydraulic_scenario')

def electrical_heating(current_A,voltage_V,flow_m3_s,density_kg_m3=786.,heat_capacity_J_kg_K=2200.,heat_fraction=1.):
    if any(not math.isfinite(v) or v<=0 for v in [voltage_V,flow_m3_s,density_kg_m3,heat_capacity_J_kg_K]):raise ValueError('Positive finite denominator inputs required')
    if not math.isfinite(current_A) or current_A<0 or not 0<=heat_fraction<=1:raise ValueError('Invalid current or heat fraction')
    return current_A*voltage_V*heat_fraction/(density_kg_m3*flow_m3_s*heat_capacity_J_kg_K)

def main():
    p=argparse.ArgumentParser();p.add_argument('--transport-json',type=Path,required=True);a=p.parse_args()
    data=json.loads(a.transport_json.read_text(encoding='utf-8'));base=data['summary'];inputs=data['inputs'];current=base['current_A'];voltage=1.85+inputs['overpotential_V']+18.5*current
    out=ROOT/'results/engineering';out.mkdir(parents=True,exist_ok=True)
    def save_csv(name,rows):
        with (out/name).open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    common=dict(length_m=inputs['length_m'],width_m=inputs['width_m'],diffusivity_m2_s=inputs['diffusivity_m2_s'])
    hyd=[hydraulic(q,height_m=h,viscosity_Pa_s=mu,**common) for q in [150.,450.,1200.] for h in [.0002,.0003,.0005] for mu in [.0002,.00035,.0006]]
    save_csv('hydraulic_scenarios.csv',hyd);standard=hydraulic(inputs['flow_rate_uL_min'],height_m=inputs['height_m'],**common)
    thermal=[dict(heat_capacity_J_kg_K=cp,electrical_heat_fraction=f,current_A=current,voltage_V=voltage,temperature_rise_K=electrical_heating(current,voltage,standard['flow_m3_s'],heat_capacity_J_kg_K=cp,heat_fraction=f),evidence_role='uncalibrated_electrical_heating_scenario') for cp in [1500.,2200.,3000.] for f in [.25,.5,1.]]
    save_csv('thermal_scenarios.csv',thermal)
    record=dict(base_hydraulics=standard,base_reviewed_current_A=current,assumed_cell_voltage_V=voltage,base_all_electrical_heat_rise_K=electrical_heating(current,voltage,standard['flow_m3_s']),hydraulic_cases=len(hyd),thermal_cases=len(thermal),transport_source_sha256=hashlib.sha256(a.transport_json.read_bytes()).hexdigest(),assumptions=['Geometry is the source 60 mm x 300 micrometre x 12 mm cell; base flow 450 uL/min.','Density 786 kg/m3 follows the source assumption; viscosity 0.00035 Pa s and Cp 2200 J/(kg K) are selected scenarios, not verified mixture properties.','Pressure estimate is fully developed parallel-plate Newtonian flow; no entrance, tubing, sidewall correction, manifold, gas evolution or pump efficiency.','Voltage = 1.85 + 0.48 + 18.5 I is the source assumed lumped voltage law; not a solved potential or measured voltage.','Temperature rise allocates a selected fraction of electrical power to sensible heat, with no heat transfer or reaction enthalpy; not a thermal PDE, measured temperature or strict safety bound.'],sources=[dict(title='MIT-hosted microfluidic mixing notes',url='https://ocw.mit.edu/courses/2-674-micro-nano-engineering-laboratory-spring-2016/f28f8672e6386d65276d2e6776ae8611_MIT2_674S16_MicrofluidcMix.pdf',supports='Definitions of Reynolds/Peclet numbers and diffusion time scaling; not fluid properties'),dict(title='COMSOL Laminar Flow interface',url='https://doc.comsol.com/6.4/doc/com.comsol.help.mfl/mfl_ug_fluidflow_single.06.04.html',supports='Shallow parallel-boundary drag and model scope; no COMSOL execution')],verification_date='2026-09-28',thermal_field_solved=False,hydraulic_CFD_solved=False,industrial_validation=False)
    (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(record,indent=2))

if __name__=='__main__':main()
