import json
G=json.load(open('guion.json')); D=json.load(open('duraciones.json'))
PROG={
 'mision':{'titulo':[0,0],'hohmann':[0,0],'ida':[0,.46],'encuentro':[.46,.50],'salto':[.50,.50],'vuelta':[.50,1],'resumen':[1,1],'datos':[1,1]},
 'encuentro':{'titulo':[0,0],'hiperbola':[0,1],'giro':[1,1],'sol':[0,1],'triangulo':[1,1],'numeros':[1,1],'energia':[1,1],'cierre':[1,1]},
 'conservacion':{'titulo':[0,0],'energia':[0,.45],'salto':[.45,.58],'momento':[.58,.75],'h':[.75,.9],'simetria':[.9,1],'voyager':[1,1]},
 'conicas':{k:[0,1] for k in ['titulo','lento','avance','curva','valores','umbral','conclusion']},
 'voyager':{k:[0,1] for k in ['titulo','datos','llegada','periapsis','puntual','maniobra','completo','cierre']},
}
DUR={
 'mision':{'titulo':2.2,'hohmann':3,'ida':6,'encuentro':2.5,'salto':3,'vuelta':6,'resumen':4},
 'encuentro':{'titulo':2.2,'hiperbola':6,'giro':3.5,'sol':6,'triangulo':4,'numeros':3,'energia':4},
 'conservacion':{'titulo':2.2,'energia':6,'salto':3,'momento':3,'h':3.5,'simetria':3.5,'voyager':4},
 'conicas':{'titulo':2.2,'lento':7,'avance':2,'curva':5,'valores':3,'umbral':3,'conclusion':3},
 'voyager':{'titulo':2.2,'datos':2,'llegada':4,'periapsis':4,'puntual':4,'maniobra':4,'completo':4,'cierre':2},
}
def build():
    cfg={}
    for name,c in G.items():
        beats=[]
        for b in c['beats']:
            if b['id'] not in DUR[name]: continue
            beats.append({'id':b['id'],'t':'','dur':DUR[name][b['id']],'prog':PROG[name][b['id']]})
        cfg[name]={'titulo':c['titulo'],'sub':c['sub'],'beats':beats}
    return cfg
if __name__=='__main__':
    c=build(); json.dump(c,open('cfg.json','w'),ensure_ascii=False)
    for n,v in c.items(): print(n,round(sum(b['dur'] for b in v['beats']),1))
