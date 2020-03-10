#, -*-, coding:, utf-8, -*-
"""
Created, on, Wed, Dec, 18, 22:50:42, 2013

@author:, Gebruiker
"""
import numpy as np


b_value=np.array([1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500,1500])
S=np.array([0.39481,0.43774,0.12879,0.31532,0.31744,0.36900,0.59490,0.35280,0.36880,0.44046,0.48088,0.17118,0.22700,0.34665,0.26000,0.25414,0.21642,0.34456,0.26625,0.20723,0.30364])
S0=1.0,
GradientOrientations=np.array([0.1639, 0.5115, 0.8435,0.1176, -0.5388, 0.8342,0.5554, 0.8278, -0.0797,-0.4804, 0.8719, 0.0948,0.9251, -0.0442, 0.3772,0.7512, -0.0273, -0.6596,0.1655, -0.0161, 0.9861,0.6129, -0.3427, 0.7120,0.6401, 0.2747, 0.7175,-0.3724, -0.3007, 0.8780,-0.3451, 0.3167, 0.8835,0.4228, 0.7872, 0.4489,0.0441, 0.9990, 0.0089,-0.1860, 0.8131, 0.5515,0.8702, 0.4606, 0.1748,-0.7239, 0.5285, 0.4434,-0.2574, -0.8032, 0.5372,0.3515, -0.8292, 0.4346,-0.7680, -0.4705, 0.4346,0.8261, -0.5384, 0.1660,0.9852, -0.0420, -0.1660])
bvecs=np.reshape(GradientOrientations, (21,3))

ydata=np.log(S/S0)

#Construct all possible monomials of a specific order
G=constructMatrixOfMonomials(GradientOrientations, order)# %computes G from section 5.1 (ISBI'10)
#Construct set of polynomial coefficients C
C=constructSetOf321Polynomials(order)# %computes C from section 5.1 (ISBI'10)
P=G*C
P=-diag(b_value)*P

def constructMatrixOfMonomials(g, order):
    G=np.empty((len(g), order+1))
    for k in range(len(g)):
        c=0
        for i in range(order):
            for j in range(order-i):
                G[k,c]=(g[k,0]**i)*(g[k,1]**j)*(g[k,2]**(order-i-j))
                c+=1
    return G
    
def construct_321_polynomials(order):
    g=UVs()
    Mprime=321
    g=g[1:Mprime,:]
    pop=np.empty((3,3,3))
    for i in range(order):
        for j in range(order-i):
            pop[i,j,order-i-j]=population(i,j,order-i-j,order)
    C=np.empty((321,3))            
    for k in range(len(g)):
        c=0
        for i in range(order):
            for j in range(order-i):
                C[k,c]=pop(i,j,order-i-j)*(g[k,1]**i)*(g[k,2]**j)*(g[k,3]**(order-i-j))
                c+=1
    return C
    
def population(i,j,k,order):
    size=3**order
    counter=0
    for z in range(size):
        c=populationBasis(z,order,np.zeros(3))
        if (c[1]==i)and(c[2]==j)and(c[3]==k):
            counter+=1
    return c

def populationBasis(i, order, c):
    if order==0:  
        ret=c
    else:
	    j=np.mod(i,3);
	    c[j]=c[j];
	    ret=populationBasis((i-j)/3,order-1,c);
    return ret
    
def UVs():
    uvs=np.array([[ 0.        ,  0.52573109,  0.85065079],
       [ 0.        , -0.52573109,  0.85065079],
       [ 0.52573109,  0.85065079,  0.        ],
       [-0.52573109,  0.85065079,  0.        ],
       [ 0.85065079,  0.        ,  0.52573109],
       [ 0.85065079,  0.        , -0.52573109],
       [ 0.        ,  0.        ,  1.        ],
       [ 0.5       , -0.309017  ,  0.809017  ],
       [ 0.5       ,  0.309017  ,  0.809017  ],
       [-0.5       , -0.309017  ,  0.809017  ],
       [-0.5       ,  0.309017  ,  0.809017  ],
       [ 0.309017  ,  0.809017  ,  0.5       ],
       [ 0.        ,  1.        ,  0.        ],
       [-0.309017  ,  0.809017  ,  0.5       ],
       [ 0.809017  ,  0.5       ,  0.309017  ],
       [-0.809017  ,  0.5       ,  0.309017  ],
       [-0.309017  , -0.809017  ,  0.5       ],
       [ 0.309017  , -0.809017  ,  0.5       ],
       [-0.809017  , -0.5       ,  0.309017  ],
       [ 0.809017  , -0.5       ,  0.309017  ],
       [ 1.        ,  0.        ,  0.        ],
       [ 0.        ,  0.27326652,  0.96193838],
       [ 0.        , -0.27326652,  0.96193838],
       [ 0.2598919 , -0.43388855,  0.86266845],
       [ 0.70204645, -0.16062205,  0.69378048],
       [ 0.2598919 ,  0.43388855,  0.86266845],
       [ 0.70204645,  0.16062205,  0.69378048],
       [-0.2598919 , -0.43388855,  0.86266845],
       [-0.70204645, -0.16062205,  0.69378048],
       [-0.2598919 ,  0.43388855,  0.86266845],
       [-0.70204645,  0.16062205,  0.69378048],
       [ 0.16062205,  0.69378048,  0.70204645],
       [ 0.43388855,  0.86266845,  0.2598919 ],
       [ 0.27326652,  0.96193838,  0.        ],
       [-0.27326652,  0.96193838,  0.        ],
       [-0.16062205,  0.69378048,  0.70204645],
       [-0.43388855,  0.86266845,  0.2598919 ],
       [ 0.69378048,  0.70204645,  0.16062205],
       [ 0.86266845,  0.2598919 ,  0.43388855],
       [-0.69378048,  0.70204645,  0.16062205],
       [-0.86266845,  0.2598919 ,  0.43388855],
       [-0.16062205, -0.69378048,  0.70204645],
       [-0.43388855, -0.86266845,  0.2598919 ],
       [ 0.16062205, -0.69378048,  0.70204645],
       [ 0.43388855, -0.86266845,  0.2598919 ],
       [-0.69378048, -0.70204645,  0.16062205],
       [-0.86266845, -0.2598919 ,  0.43388855],
       [ 0.69378048, -0.70204645,  0.16062205],
       [ 0.86266845, -0.2598919 ,  0.43388855],
       [ 0.96193838,  0.        ,  0.27326652],
       [ 0.96193838,  0.        , -0.27326652],
       [ 0.26286557, -0.16245985,  0.95105654],
       [ 0.52573109,  0.        ,  0.85065079],
       [ 0.26286557,  0.16245985,  0.95105654],
       [-0.26286557, -0.16245985,  0.95105654],
       [-0.52573109,  0.        ,  0.85065079],
       [-0.26286557,  0.16245985,  0.95105654],
       [ 0.16245985,  0.95105654,  0.26286557],
       [-0.16245985,  0.95105654,  0.26286557],
       [ 0.        ,  0.85065079,  0.52573109],
       [ 0.58778524,  0.68819094,  0.42532539],
       [ 0.68819094,  0.42532539,  0.58778524],
       [ 0.42532539,  0.58778524,  0.68819094],
       [-0.58778524,  0.68819094,  0.42532539],
       [-0.68819094,  0.42532539,  0.58778524],
       [-0.42532539,  0.58778524,  0.68819094],
       [-0.16245985, -0.95105654,  0.26286557],
       [ 0.16245985, -0.95105654,  0.26286557],
       [ 0.        , -0.85065079,  0.52573109],
       [-0.58778524, -0.68819094,  0.42532539],
       [-0.68819094, -0.42532539,  0.58778524],
       [-0.42532539, -0.58778524,  0.68819094],
       [ 0.58778524, -0.68819094,  0.42532539],
       [ 0.68819094, -0.42532539,  0.58778524],
       [ 0.42532539, -0.58778524,  0.68819094],
       [ 0.95105654,  0.26286557,  0.16245985],
       [ 0.95105654,  0.26286557, -0.16245985],
       [ 0.85065079,  0.52573109,  0.        ],
       [-0.95105654,  0.26286557, -0.16245985],
       [-0.95105654,  0.26286557,  0.16245985],
       [-0.85065079,  0.52573109,  0.        ],
       [ 0.        ,  0.40335536,  0.91504341],
       [ 0.        ,  0.13795224,  0.99043888],
       [ 0.        , -0.40335536,  0.91504341],
       [ 0.        , -0.13795224,  0.99043888],
       [ 0.13120037, -0.48444164,  0.86492932],
       [ 0.38361374, -0.37503856,  0.84391147],
       [ 0.78384304, -0.08108629,  0.61564201],
       [ 0.60682517, -0.23708633,  0.75865233],
       [ 0.13120037,  0.48444164,  0.86492932],
       [ 0.38361374,  0.37503856,  0.84391147],
       [ 0.78384304,  0.08108629,  0.61564201],
       [ 0.60682517,  0.23708633,  0.75865233],
       [-0.13120037, -0.48444164,  0.86492932],
       [-0.38361374, -0.37503856,  0.84391147],
       [-0.78384304, -0.08108629,  0.61564201],
       [-0.60682517, -0.23708633,  0.75865233],
       [-0.13120037,  0.48444164,  0.86492932],
       [-0.38361374,  0.37503856,  0.84391147],
       [-0.78384304,  0.08108629,  0.61564201],
       [-0.60682517,  0.23708633,  0.75865233],
       [ 0.08108629,  0.61564201,  0.78384304],
       [ 0.23708633,  0.75865233,  0.60682517],
       [ 0.48444164,  0.86492932,  0.13120037],
       [ 0.37503856,  0.84391147,  0.38361374],
       [ 0.40335536,  0.91504341,  0.        ],
       [ 0.13795224,  0.99043888,  0.        ],
       [-0.40335536,  0.91504341,  0.        ],
       [-0.13795224,  0.99043888,  0.        ],
       [-0.08108629,  0.61564201,  0.78384304],
       [-0.23708633,  0.75865233,  0.60682517],
       [-0.48444164,  0.86492932,  0.13120037],
       [-0.37503856,  0.84391147,  0.38361374],
       [ 0.61564201,  0.78384304,  0.08108629],
       [ 0.75865233,  0.60682517,  0.23708633],
       [ 0.86492932,  0.13120037,  0.48444164],
       [ 0.84391147,  0.38361374,  0.37503856],
       [-0.61564201,  0.78384304,  0.08108629],
       [-0.75865233,  0.60682517,  0.23708633],
       [-0.86492932,  0.13120037,  0.48444164],
       [-0.84391147,  0.38361374,  0.37503856],
       [-0.08108629, -0.61564201,  0.78384304],
       [-0.23708633, -0.75865233,  0.60682517],
       [-0.48444164, -0.86492932,  0.13120037],
       [-0.37503856, -0.84391147,  0.38361374],
       [ 0.08108629, -0.61564201,  0.78384304],
       [ 0.23708633, -0.75865233,  0.60682517],
       [ 0.48444164, -0.86492932,  0.13120037],
       [ 0.37503856, -0.84391147,  0.38361374],
       [-0.61564201, -0.78384304,  0.08108629],
       [-0.75865233, -0.60682517,  0.23708633],
       [-0.86492932, -0.13120037,  0.48444164],
       [-0.84391147, -0.38361374,  0.37503856],
       [ 0.61564201, -0.78384304,  0.08108629],
       [ 0.75865233, -0.60682517,  0.23708633],
       [ 0.86492932, -0.13120037,  0.48444164],
       [ 0.84391147, -0.38361374,  0.37503856],
       [ 0.91504341,  0.        ,  0.40335536],
       [ 0.99043888,  0.        ,  0.13795224],
       [ 0.91504341,  0.        , -0.40335536],
       [ 0.99043888,  0.        , -0.13795224],
       [ 0.13307109, -0.08224247,  0.98768836],
       [ 0.3861874 , -0.23867694,  0.89100653],
       [ 0.5192585 , -0.15643448,  0.84017789],
       [ 0.5192585 ,  0.15643448,  0.84017789],
       [ 0.13307109,  0.08224247,  0.98768836],
       [ 0.3861874 ,  0.23867694,  0.89100653],
       [-0.13307109, -0.08224247,  0.98768836],
       [-0.3861874 , -0.23867694,  0.89100653],
       [-0.5192585 , -0.15643448,  0.84017789],
       [-0.5192585 ,  0.15643448,  0.84017789],
       [-0.13307109,  0.08224247,  0.98768836],
       [-0.3861874 ,  0.23867694,  0.89100653],
       [ 0.23867694,  0.89100653,  0.3861874 ],
       [ 0.08224247,  0.98768836,  0.13307109],
       [-0.08224247,  0.98768836,  0.13307109],
       [-0.23867694,  0.89100653,  0.3861874 ],
       [ 0.15643448,  0.84017789,  0.5192585 ],
       [-0.15643448,  0.84017789,  0.5192585 ],
       [ 0.45399049,  0.7579354 ,  0.46842986],
       [ 0.70710677,  0.60150099,  0.37174803],
       [ 0.7579354 ,  0.46842986,  0.45399049],
       [ 0.60150099,  0.37174803,  0.70710677],
       [ 0.37174803,  0.70710677,  0.60150099],
       [ 0.46842986,  0.45399049,  0.7579354 ],
       [-0.45399049,  0.7579354 ,  0.46842986],
       [-0.70710677,  0.60150099,  0.37174803],
       [-0.7579354 ,  0.46842986,  0.45399049],
       [-0.60150099,  0.37174803,  0.70710677],
       [-0.37174803,  0.70710677,  0.60150099],
       [-0.46842986,  0.45399049,  0.7579354 ],
       [-0.23867694, -0.89100653,  0.3861874 ],
       [-0.08224247, -0.98768836,  0.13307109],
       [ 0.08224247, -0.98768836,  0.13307109],
       [ 0.23867694, -0.89100653,  0.3861874 ],
       [-0.15643448, -0.84017789,  0.5192585 ],
       [ 0.15643448, -0.84017789,  0.5192585 ],
       [-0.45399049, -0.7579354 ,  0.46842986],
       [-0.70710677, -0.60150099,  0.37174803],
       [-0.7579354 , -0.46842986,  0.45399049],
       [-0.60150099, -0.37174803,  0.70710677],
       [-0.37174803, -0.70710677,  0.60150099],
       [-0.46842986, -0.45399049,  0.7579354 ],
       [ 0.45399049, -0.7579354 ,  0.46842986],
       [ 0.70710677, -0.60150099,  0.37174803],
       [ 0.7579354 , -0.46842986,  0.45399049],
       [ 0.60150099, -0.37174803,  0.70710677],
       [ 0.37174803, -0.70710677,  0.60150099],
       [ 0.46842986, -0.45399049,  0.7579354 ],
       [ 0.89100653,  0.3861874 ,  0.23867694],
       [ 0.98768836,  0.13307109,  0.08224247],
       [ 0.98768836,  0.13307109, -0.08224247],
       [ 0.89100653,  0.3861874 , -0.23867694],
       [ 0.84017789,  0.5192585 ,  0.15643448],
       [ 0.84017789,  0.5192585 , -0.15643448],
       [-0.89100653,  0.3861874 , -0.23867694],
       [-0.98768836,  0.13307109, -0.08224247],
       [-0.98768836,  0.13307109,  0.08224247],
       [-0.89100653,  0.3861874 ,  0.23867694],
       [-0.84017789,  0.5192585 , -0.15643448],
       [-0.84017789,  0.5192585 ,  0.15643448],
       [ 0.13165537, -0.3582288 ,  0.9243046 ],
       [ 0.26408276, -0.30125889,  0.91624421],
       [ 0.13279247, -0.22011703,  0.96639258],
       [ 0.71128172,  0.        ,  0.70290703],
       [ 0.62023956,  0.08114185,  0.78020436],
       [ 0.62023956, -0.08114185,  0.78020436],
       [ 0.13165537,  0.3582288 ,  0.9243046 ],
       [ 0.26408276,  0.30125889,  0.91624421],
       [ 0.13279247,  0.22011703,  0.96639258],
       [ 0.39960706, -0.08232358,  0.91298246],
       [ 0.39960706,  0.08232358,  0.91298246],
       [ 0.26640469,  0.        ,  0.96386129],
       [-0.13165537, -0.3582288 ,  0.9243046 ],
       [-0.26408276, -0.30125889,  0.91624421],
       [-0.13279247, -0.22011703,  0.96639258],
       [-0.71128172,  0.        ,  0.70290703],
       [-0.62023956,  0.08114185,  0.78020436],
       [-0.62023956, -0.08114185,  0.78020436],
       [-0.13165537,  0.3582288 ,  0.9243046 ],
       [-0.26408276,  0.30125889,  0.91624421],
       [-0.13279247,  0.22011703,  0.96639258],
       [-0.39960706, -0.08232358,  0.91298246],
       [-0.39960706,  0.08232358,  0.91298246],
       [-0.26640469,  0.        ,  0.96386129],
       [ 0.3582288 ,  0.9243046 ,  0.13165537],
       [ 0.22011703,  0.96639258,  0.13279247],
       [ 0.30125889,  0.91624421,  0.26408276],
       [-0.3582288 ,  0.9243046 ,  0.13165537],
       [-0.30125889,  0.91624421,  0.26408276],
       [-0.22011703,  0.96639258,  0.13279247],
       [ 0.        ,  0.70290703,  0.71128172],
       [-0.08114185,  0.78020436,  0.62023956],
       [ 0.08114185,  0.78020436,  0.62023956],
       [ 0.        ,  0.96386129,  0.26640469],
       [-0.08232358,  0.91298246,  0.39960706],
       [ 0.08232358,  0.91298246,  0.39960706],
       [ 0.57125163,  0.79264921,  0.21302287],
       [ 0.64741188,  0.70230985,  0.29600459],
       [ 0.51612163,  0.78345168,  0.34615302],
       [ 0.79264921,  0.21302287,  0.57125163],
       [ 0.70230985,  0.29600459,  0.64741188],
       [ 0.78345168,  0.34615302,  0.51612163],
       [ 0.21302287,  0.57125163,  0.79264921],
       [ 0.34615302,  0.51612163,  0.78345168],
       [ 0.29600459,  0.64741188,  0.70230985],
       [ 0.64657778,  0.56425422,  0.51337546],
       [ 0.56425422,  0.51337546,  0.64657778],
       [ 0.51337546,  0.64657778,  0.56425422],
       [-0.57125163,  0.79264921,  0.21302287],
       [-0.64741188,  0.70230985,  0.29600459],
       [-0.51612163,  0.78345168,  0.34615302],
       [-0.79264921,  0.21302287,  0.57125163],
       [-0.70230985,  0.29600459,  0.64741188],
       [-0.78345168,  0.34615302,  0.51612163],
       [-0.21302287,  0.57125163,  0.79264921],
       [-0.34615302,  0.51612163,  0.78345168],
       [-0.29600459,  0.64741188,  0.70230985],
       [-0.64657778,  0.56425422,  0.51337546],
       [-0.56425422,  0.51337546,  0.64657778],
       [-0.51337546,  0.64657778,  0.56425422],
       [-0.3582288 , -0.9243046 ,  0.13165537],
       [-0.22011703, -0.96639258,  0.13279247],
       [-0.30125889, -0.91624421,  0.26408276],
       [ 0.3582288 , -0.9243046 ,  0.13165537],
       [ 0.30125889, -0.91624421,  0.26408276],
       [ 0.22011703, -0.96639258,  0.13279247],
       [ 0.        , -0.70290703,  0.71128172],
       [ 0.08114185, -0.78020436,  0.62023956],
       [-0.08114185, -0.78020436,  0.62023956],
       [ 0.        , -0.96386129,  0.26640469],
       [ 0.08232358, -0.91298246,  0.39960706],
       [-0.08232358, -0.91298246,  0.39960706],
       [-0.57125163, -0.79264921,  0.21302287],
       [-0.64741188, -0.70230985,  0.29600459],
       [-0.51612163, -0.78345168,  0.34615302],
       [-0.79264921, -0.21302287,  0.57125163],
       [-0.70230985, -0.29600459,  0.64741188],
       [-0.78345168, -0.34615302,  0.51612163],
       [-0.21302287, -0.57125163,  0.79264921],
       [-0.34615302, -0.51612163,  0.78345168],
       [-0.29600459, -0.64741188,  0.70230985],
       [-0.64657778, -0.56425422,  0.51337546],
       [-0.56425422, -0.51337546,  0.64657778],
       [-0.51337546, -0.64657778,  0.56425422],
       [ 0.57125163, -0.79264921,  0.21302287],
       [ 0.64741188, -0.70230985,  0.29600459],
       [ 0.51612163, -0.78345168,  0.34615302],
       [ 0.79264921, -0.21302287,  0.57125163],
       [ 0.70230985, -0.29600459,  0.64741188],
       [ 0.78345168, -0.34615302,  0.51612163],
       [ 0.21302287, -0.57125163,  0.79264921],
       [ 0.34615302, -0.51612163,  0.78345168],
       [ 0.29600459, -0.64741188,  0.70230985],
       [ 0.64657778, -0.56425422,  0.51337546],
       [ 0.56425422, -0.51337546,  0.64657778],
       [ 0.51337546, -0.64657778,  0.56425422],
       [ 0.9243046 ,  0.13165537,  0.3582288 ],
       [ 0.96639258,  0.13279247,  0.22011703],
       [ 0.91624421,  0.26408276,  0.30125889],
       [ 0.9243046 ,  0.13165537, -0.3582288 ],
       [ 0.91624421,  0.26408276, -0.30125889],
       [ 0.96639258,  0.13279247, -0.22011703],
       [ 0.70290703,  0.71128172,  0.        ],
       [ 0.78020436,  0.62023956, -0.08114185],
       [ 0.78020436,  0.62023956,  0.08114185],
       [ 0.96386129,  0.26640469,  0.        ],
       [ 0.91298246,  0.39960706, -0.08232358],
       [ 0.91298246,  0.39960706,  0.08232358],
       [-0.9243046 ,  0.13165537, -0.3582288 ],
       [-0.96639258,  0.13279247, -0.22011703],
       [-0.91624421,  0.26408276, -0.30125889],
       [-0.9243046 ,  0.13165537,  0.3582288 ],
       [-0.91624421,  0.26408276,  0.30125889],
       [-0.96639258,  0.13279247,  0.22011703],
       [-0.70290703,  0.71128172,  0.        ],
       [-0.78020436,  0.62023956,  0.08114185],
       [-0.78020436,  0.62023956, -0.08114185],
       [-0.96386129,  0.26640469,  0.        ],
       [-0.91298246,  0.39960706,  0.08232358],
       [-0.91298246,  0.39960706, -0.08232358]])
    return uvs