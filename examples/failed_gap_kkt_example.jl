using PolyPlanning
using SparseArrays

maps=[
[ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999964633505863, 0.9999964633505863, 0.0026595650621026607, 1.0, 0.0026595650621026607, -1.0], 4, 2), 
[2.6241928467766584, -0.3, -2.375789469976273, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.3], [2.375, 0.3]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 
2], [0.9999964633505863, -0.9999964633505863, 1.0, -0.0026595650621026607, -1.0, -0.0026595650621026607], 4, 2), [5.0, -2.3757894699762736, -0.3, 2.6241928467766584], [[2.375, -0.3], [2.625, -0.3], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999964474556374, 0.9999964474556374, 0.0026655348627465463, 1.0, 0.0026655348627465463, -1.0], 4, 2), 
[2.6241629558505113, -0.3105263157894737, -2.3758192814276757, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.3105263157894737], 
[2.375, 0.3105263157894737]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 
2, 2], [0.9999964474556374, -0.9999964474556374, 1.0, -0.0026655348627465463, -1.0, -0.0026655348627465463], 4, 2), [5.0, -2.375819281427676, -0.3105263157894737, 2.6241629558505113], [[2.375, -0.3105263157894737], 
[2.625, -0.3105263157894737], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999964314532913, 0.9999964314532913, 0.00267153152384882, 1.0, 0.00267153152384882, -1.0], 4, 2), [2.624132930338812, -0.32105263157894737, -2.3758492269276443, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.32105263157894737], [2.375, 0.32105263157894737]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999964314532913, -0.9999964314532913, 1.0, -0.00267153152384882, -1.0, -0.00267153152384882], 4, 2), [5.0, -2.3758492269276448, -0.32105263157894737, 2.624132930338812], [[2.375, -0.32105263157894737], [2.625, -0.32105263157894737], [2.6375, -5.0], 
[2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999964153425781, 0.9999964153425781, 0.0026775552271010787, 1.0, 0.0026775552271010787, -1.0], 4, 2), 
[2.6241027693305443, -0.33157894736842103, 
-2.3758793073823457, 5.0], [[2.3625, 5.0], 
[2.6375, 5.0], [2.625, 0.33157894736842103], [2.375, 0.33157894736842103]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 
2, 2, 2], [0.9999964153425781, -0.9999964153425781, 1.0, -0.0026775552271010787, -1.0, -0.0026775552271010787], 4, 2), [5.0, -2.375879307382346, -0.33157894736842103, 2.6241027693305448], [[2.375, -0.33157894736842103], [2.625, -0.33157894736842103], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.999996399122517, 
0.999996399122517, 0.0026836061558373012, 1.0, 0.0026836061558373012, -1.0], 4, 2), [2.6240724719064525, -0.34210526315789475, -2.375909523706133, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.34210526315789475], [2.375, 0.34210526315789475]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 
2, 2], [0.999996399122517, -0.999996399122517, 1.0, -0.0026836061558373012, -1.0, -0.0026836061558373012], 4, 2), [5.0, -2.3759095237061327, -0.34210526315789475, 2.624072471906452], [[2.375, -0.34210526315789475], 
[2.625, -0.34210526315789475], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999963827921166, 0.9999963827921166, 0.002689684495052447, 
1.0, 0.002689684495052447, -1.0], 4, 2), [2.6240420371389455, -0.3526315789473684, -2.3759398768216373, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.3526315789473684], [2.375, 0.3526315789473684]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999963827921166, -0.9999963827921166, 1.0, -0.002689684495052447, -1.0, -0.002689684495052447], 4, 2), [5.0, -2.3759398768216373, -0.3526315789473684, 2.6240420371389455], [[2.375, -0.3526315789473684], [2.625, -0.3526315789473684], [2.6375, -5.0], 
[2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999963663503731, 0.9999963663503731, 0.0026957904314213054, 1.0, 0.0026957904314213054, -1.0], 4, 2), 
[2.624011464092003, -0.3631578947368421, -2.3759703676598627, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.3631578947368421], [2.375, 0.3631578947368421]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999963663503731, -0.9999963663503731, 1.0, -0.0026957904314213054, -1.0, -0.0026957904314213054], 4, 2), [5.0, -2.375970367659863, -0.3631578947368421, 2.624011464092003], [[2.375, -0.3631578947368421], [2.625, -0.3631578947368421], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999963497962729, 0.9999963497962729, 0.002701924153317613, 
1.0, 0.002701924153317613, -1.0], 4, 2), [2.6239807518210823, -0.3736842105263158, -2.3760009971602827, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.3736842105263158], [2.375, 0.3736842105263158]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999963497962729, -0.9999963497962729, 1.0, -0.002701924153317613, -1.0, -0.002701924153317613], 4, 2), [5.0, -2.3760009971602827, -0.3736842105263158, 2.6239807518210823], [[2.375, -0.3736842105263158], [2.625, -0.3736842105263158], [2.6375, -5.0], 
[2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999963331287892, 0.9999963331287892, 0.0027080858508334187, 1.0, 0.0027080858508334187, -1.0], 4, 2), 
[2.6239498993730144, -0.38421052631578945, 
-2.376031766270932, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.38421052631578945], [2.375, 0.38421052631578945]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999963331287892, -0.9999963331287892, 1.0, -0.0027080858508334187, -1.0, 
-0.0027080858508334187], 4, 2), [5.0, -2.3760317662709314, -0.38421052631578945, 2.6239498993730144], [[2.375, -0.38421052631578945], [2.625, -0.38421052631578945], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999963163468847, 0.9999963163468847, 0.0027142757157987253, 1.0, 0.0027142757157987253, -1.0], 4, 2), 
[2.623918905785915, -0.39473684210526316, -2.3760626759485084, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.39473684210526316], [2.375, 0.39473684210526316]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999963163468847, -0.9999963163468847, 1.0, -0.0027142757157987253, -1.0, 
-0.0027142757157987253], 4, 2), [5.0, -2.3760626759485084, -0.39473684210526316, 2.623918905785915], [[2.375, -0.39473684210526316], [2.625, -0.39473684210526316], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999962994495093, 0.9999962994495093, 0.0027204939418013956, 1.0, 0.0027204939418013956, -1.0], 4, 2), 
[2.623887770089074, -0.4052631578947368, -2.376093727158472, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.4052631578947368], [2.375, 0.4052631578947368]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999962994495093, -0.9999962994495093, 1.0, -0.0027204939418013956, -1.0, -0.0027204939418013956], 4, 2), [5.0, -2.3760937271584726, -0.4052631578947368, 2.623887770089074], [[2.375, -0.4052631578947368], [2.625, -0.4052631578947368], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999962824356013, 0.9999962824356013, 0.0027267407242073326, 1.0, 0.0027267407242073326, -1.0], 4, 2), 
[2.6238564913028624, -0.41578947368421054, 
-2.376124920875145, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.41578947368421054], [2.375, 0.41578947368421054]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999962824356013, -0.9999962824356013, 1.0, -0.0027267407242073326, -1.0, 
-0.0027267407242073326], 4, 2), [5.0, -2.3761249208751445, -0.41578947368421054, 2.623856491302862], [[2.375, -0.41578947368421054], [2.625, -0.41578947368421054], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999962653040869, 0.9999962653040869, 0.0027330162601809436, 1.0, 0.0027330162601809436, -1.0], 4, 2), 
[2.6238250684386246, -0.4263157894736842, -2.37615625808181, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.4263157894736842], [2.375, 0.4263157894736842]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999962653040869, -0.9999962653040869, 1.0, -0.0027330162601809436, -1.0, -0.0027330162601809436], 4, 2), [5.0, -2.37615625808181, -0.4263157894736842, 2.6238250684386246], [[2.375, -0.4263157894736842], [2.625, -0.4263157894736842], [2.6375, -5.0], 
[2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999962480538794, 0.9999962480538794, 0.002739320748705879, 
1.0, 0.002739320748705879, -1.0], 4, 2), [2.6237935004985777, -0.4368421052631579, -2.3761877397708195, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.4368421052631579], [2.375, 0.4368421052631579]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999962480538794, -0.9999962480538794, 1.0, -0.002739320748705879, -1.0, -0.002739320748705879], 4, 2), [5.0, -2.376187739770819, -0.4368421052631579, 2.6237935004985777], [[2.375, -0.4368421052631579], [2.625, -0.4368421052631579], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999962306838799, 0.9999962306838799, 0.0027456543906060676, 1.0, 0.0027456543906060676, -1.0], 4, 2), 
[2.6237617864757032, -0.4473684210526316, -2.3762193669436966, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.4473684210526316], 
[2.375, 0.4473684210526316]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 
2, 2], [0.9999962306838799, -0.9999962306838799, 1.0, -0.0027456543906060676, -1.0, -0.0027456543906060676], 4, 2), [5.0, -2.376219366943696, -0.4473684210526316, 2.623761786475703], [[2.375, -0.4473684210526316], [2.625, -0.4473684210526316], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999962131929765, 0.9999962131929765, 0.002752017388567037, 
1.0, 0.002752017388567037, -1.0], 4, 2), [2.6237299253536404, -0.45789473684210524, -2.376251140611242, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.45789473684210524], [2.375, 0.45789473684210524]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 
2, 2], [0.9999962131929765, -0.9999962131929765, 1.0, -0.002752017388567037, -1.0, -0.002752017388567037], 4, 2), [5.0, -2.376251140611242, -0.45789473684210524, 2.6237299253536404], [[2.375, -0.45789473684210524], 
[2.625, -0.45789473684210524], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.9999961955800449, 0.9999961955800449, 0.0027584099471575385, 1.0, 0.0027584099471575385, -1.0], 4, 2), 
[2.623697916106581, -0.46842105263157896, -2.3762830617936435, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.46842105263157896], [2.375, 0.46842105263157896]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.9999961955800449, -0.9999961955800449, 1.0, -0.0027584099471575385, -1.0, 
-0.0027584099471575385], 4, 2), [5.0, -2.3762830617936435, -0.46842105263157896, 2.623697916106581], [[2.375, -0.46842105263157896], [2.625, -0.46842105263157896], [2.6375, -5.0], [2.3625, -5.0]])],

 [ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.999996177843947, 
0.999996177843947, 0.0027648322728514646, 1.0, 0.0027648322728514646, -1.0], 4, 2), [2.623665757699153, -0.4789473684210526, -2.376315131520582, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.4789473684210526], [2.375, 0.4789473684210526]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.999996177843947, -0.999996177843947, 
1.0, -0.0027648322728514646, -1.0, -0.0027648322728514646], 4, 2), [5.0, -2.376315131520582, -0.4789473684210526, 2.623665757699153], [[2.375, -0.4789473684210526], [2.625, -0.4789473684210526], [2.6375, -5.0], [2.3625, -5.0]])],

[ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.999996159983532, 
0.999996159983532, 0.0027712845740500845, 1.0, 0.0027712845740500845, -1.0], 4, 2), [2.623633449086315, -0.48947368421052634, -2.3763473508313444, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.48947368421052634], [2.375, 0.48947368421052634]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 
2, 2], [0.999996159983532, -0.999996159983532, 1.0, -0.0027712845740500845, -1.0, -0.0027712845740500845], 4, 2), [5.0, -2.3763473508313444, -0.48947368421052634, 2.6236334490863156], [[2.375, -0.48947368421052634], [2.625, -0.48947368421052634], [2.6375, -5.0], [2.3625, -5.0]])],

[ConvexPolygon2D(sparse([1, 3, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [-0.999996141997635, 
0.999996141997635, 0.002777767061104581, 1.0, 0.002777767061104581, -1.0], 4, 2), [2.6236009892132395, -0.5, -2.3763797207749353, 5.0], [[2.3625, 5.0], [2.6375, 5.0], [2.625, 0.5], [2.375, 0.5]]) ConvexPolygon2D(sparse([2, 4, 1, 2, 3, 4], [1, 1, 2, 2, 2, 2], [0.999996141997635, -0.999996141997635, 1.0, -0.002777767061104581, -1.0, -0.002777767061104581], 4, 2), [5.0, -2.3763797207749358, -0.5, 2.6236009892132395], [[2.375, -0.5], [2.625, -0.5], [2.6375, -5.0], [2.3625, -5.0]])]
]


x0s=[
[4.514598145403474, -0.07041599465061932, 
2.8687742761945243, 0.0, 0.0, 0.0],
[5.656800299565696, -0.06436869027706393, 
2.9061925123451395, 0.0, 0.0, 0.0],
[5.987134257419568, -0.6575002532477763, -2.4673806078495235, 0.0, 0.0, 0.0],
[5.443419548101982, 0.8098546343455364, -1.0402996362311043, 0.0, 0.0, 0.0],
[4.687021710221183, -0.8763768879626046, -1.4303323086230353, 0.0, 0.0, 0.0],
[4.070565758246018, 0.3589013178289575, -2.6204344006100873, 0.0, 0.0, 0.0],
[4.563586329982898, 0.8786980024118285, -2.948811094963895, 0.0, 0.0, 0.0],
[4.342149278480869, 0.17308174709596758, 1.5781080835468622, 0.0, 0.0, 0.0],
[5.4543245551091895, 0.29214846637885694, 
0.5255385271905229, 0.0, 0.0, 0.0],
[4.451543962649737, -0.9011912853135453, -0.995072782693657, 0.0, 0.0, 0.0],
[4.669766066619071, 0.8128344780324899, -2.0812776468011744, 0.0, 0.0, 0.0],
[4.103624880944725, 0.8070808836203995, -0.6665505789082422, 0.0, 0.0, 0.0],
[5.68598944830744, -0.22865061436223555, 1.9237958830483288, 0.0, 0.0, 0.0],
[5.777545986920867, -0.14423684348856947, 
1.1945894893096414, 0.0, 0.0, 0.0],
[5.0595721838145105, 0.5782312336292672, -1.898572816485293, 0.0, 0.0, 0.0],
[5.186097430635677, -0.07142072521878351, 
1.0708017987715648, 0.0, 0.0, 0.0],
[4.61048790416723, -0.9075926792953479, 1.7445944723856934, 0.0, 0.0, 0.0],
[4.3137748043815485, -0.508384908060481, 1.0317628706457622, 0.0, 0.0, 0.0],
[5.81145897409044, -0.8696440270261034, 1.316255351254914, 0.0, 0.0, 0.0],
[4.3194870289120715, 0.6977051162626944, -2.39059118995388, 0.0, 0.0, 0.0]]
ego_rect = [PolyPlanning.ConvexPolygon2D(sparse([2, 4, 1, 3], [1, 1, 2, 2], [1.0, -1.0, 1.0, -1.0], 4, 2), [0.25, 1.0, 0.25, 1.0], [[-1.0, 0.25], [1.0, 0.25], [1.0, -0.25], [-1.0, -0.25]]) ]
Rf = 1e-3 * PolyPlanning.I(3);
Rf[3, 3] = Rf[3, 3] / 100.0;


obs_polys=maps[3]

kkt_prob = PolyPlanning.setup_direct_kkt(
    ego_rect,
    obs_polys;
    T=20,
    dt=0.2,
    Rf,
    Qf=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π
);

nonsmooth_prob = PolyPlanning.setup_nonsmooth(
    ego_rect,
    obs_polys;
    T=20,
    dt=.2,
    R_cost=Rf,
    Q_cost=2e-3 * PolyPlanning.I(2),
    u1_max=10.0,
    u2_max=10.0,
    u3_max=π,
	n_sd_slots=4
)


# x0s[11] failed
x0 = x0s[11]
sol_kkt = PolyPlanning.solve_prob_direct_kkt(kkt_prob, x0; is_displaying=false, sleep_duration=0.2)















