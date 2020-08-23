from sklearn.decomposition import PCA
from sklearn.manifold import Isomap, MDS, TSNE, LocallyLinearEmbedding as LLE
from umap import UMAP 
import numpy as np
import random
import math
import datetime
import json
import multiprocessing

max_D = 8;
max_N = 32;

#max_D = 2;
#max_N = 3;

ms = datetime.timedelta(microseconds=1000)

def benchmark(DR, key):
    #warmup
    '''print("start warmup")
    for _ in range(100):
        X = np.array([[random.random() for _ in range(50)] for _ in range(50)])
        dr = DR(n_components=2)
        if key == "TSNE" or key == "UMAP":
            dr.n_iter = 300
            if key == "TSNE":
                dr.method = "exact"
        dr.fit_transform(X=X)
    print("end warmup")'''
    finish_row = True

    B = []
    for D in range(max_D):
        d = math.floor(2 * 5**((D+1)/2))
        finish_row = True
        Brow =[]
        for N in range(max_N):
            if not finish_row:
                Brow.append(0)
                break;
            n = math.floor(16 + 2 ** (N/2 + 1))
            X = np.array([[random.random() for _ in range(d)] for _ in range(n)])
            dur = []
            for i in range(5):
                #m = []
                start = datetime.datetime.now()
                #m.append(start)
                dr = DR(n_components=2)
                if key == "TSNE" or key == "UMAP":
                    dr.n_iter = 350
                    if key == "TSNE":
                        dr.method = "exact"
                #m.append(datetime.datetime.now())
                t = True
                
                if not (key == "TSNE" or key == "MDS"):
                    dr.fit(X=X)
                mid = datetime.datetime.now()
                def f():
                    if key == "TSNE" or key == "MDS":
                        dr.fit_transform(X=X)
                    else:
                        dr.transform(X=X)
                p = multiprocessing.Process(target=f, name="DR")
                p.start()
                p.join(10)
                if p.is_alive():
                    t = False
                    p.terminate()
                    p.join()
                ende = datetime.datetime.now()
                #m.append(datetime.datetime.now())
                m = [0, (mid-start) / ms, (ende-start) / ms]
                dur.append(m)
                
                if (m[2]-m[0]) > 10000:
                    finish_row = False
                del dr  

            print(key, d, n, list(map(lambda d: d[2] - d[0], dur)), t) 
            del X 
            Brow.append({
                "d": d, 
                "n": n, 
                "dur": dur,
                "terminated": t,
            })
        B.append(Brow)
        print("\n")
    return B
        

""" X = np.array([[random.random() for _ in range(50)] for _ in range(50)])
#print(X)
start = datetime.datetime.now()
dr = TSNE(n_components=2)
dr.n_iter = 300
mid = datetime.datetime.now()
dr.fit_transform(X=X)
print((datetime.datetime.now() - start).microseconds / 1000, (mid - start).microseconds / 1000)
 """
DRS = [("PCA", PCA), ("LLE", LLE), ("TSNE", TSNE), ("MDS", MDS), ("ISOMAP", Isomap), ("UMAP", UMAP)];
DRS = [("TSNE", TSNE), ("MDS", MDS)];
data = {key: benchmark(dr, key) for key, dr in DRS}

with open('eval_sklearn_2.json', 'w') as outfile:
    json.dump(data, outfile)

