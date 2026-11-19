import json
def load(start_json,fit_json):
    P=json.load(open(start_json)); P.pop('_x',None); r=json.load(open(fit_json)); Q=dict(P); Q.update(dict(zip(r['keys'],r['x']))); return P,Q
