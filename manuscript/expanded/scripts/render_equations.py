"""Render the supplied display equations with existing Matplotlib mathtext.

This bounded fallback supports these specific expressions; it is not a general
LaTeX engine. Equation source is retained and each final PDF is visually checked.
"""
import json
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.mathtext import MathTextParser
from PIL import Image

BASE=Path(__file__).resolve().parents[1]/'assets'/'equations'
plt.rcParams.update({'mathtext.fontset':'stix','font.family':'STIXGeneral','svg.fonttype':'path'})

def main():
    expressions=json.loads((BASE/'equations.json').read_text('utf-8'))
    records={}
    parser=MathTextParser('path')
    for key,tex in expressions.items():
        clean=tex.replace(r'\begin{aligned}','').replace(r'\end{aligned}','').replace('&','')
        lines=[' '.join(s.split()) for s in re.split(r'\\\\',clean) if s.strip()]
        fig=plt.figure(figsize=(7,.32*len(lines)))
        for i,line in enumerate(lines):
            assert '\n' not in line
            parser.parse('$'+line+'$',dpi=72)
            fig.text(.5,(len(lines)-i-.5)/len(lines),'$'+line+'$',fontsize=12,ha='center',va='center')
        for extension in ['png','svg']:
            fig.savefig(BASE/(key+'.'+extension),dpi=600,bbox_inches='tight',pad_inches=.06,transparent=True)
        plt.close(fig)
        with Image.open(BASE/(key+'.png')) as im:width=im.width/600
        records[key]={'width_inches':min(6.5,width),'tex':tex,'mathtext_parser_validated':True,'renderer':'Matplotlib '+matplotlib.__version__+' mathtext/STIX; display image'}
    (BASE/'equations_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'equations':len(records),'renderer':'Matplotlib mathtext'}))

if __name__=='__main__':main()
