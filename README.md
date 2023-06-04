# Demonstration code of the article:

Cascaded Sobol' Sampling, Loïs Paulin, David Coeurjolly, Jean-Claude Iehl, Nicolas Bonneel, Alexander Keller, Victor Ostromoukhov, ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia), 40(6), pp. 274:1–274:13, December 2021.

![](https://projet.liris.cnrs.fr/cascaded/teaser.png)


Building
========

To build project:

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

Running
=======

Generate a point set of 16 points in 3 dimension in file `owenCascadedSobol_3D_16pts.dat` using `cascaded_sobol_init_tab.dat` as sobol generator matrices and with owen scrambling

    ./cascadedSobol -d 3 -n 16 -i ../data/cascaded_sobol_init_tab.dat -o owenCascadedSobol_3D_16pts.dat -p

License
=======

```
 MIT License
 
 Copyright (c) 2021 CNRS
 
 Loïs Paulin, David Coeurjolly, Jean-Claude Iehl, Nicolas Bonneel, Alexander Keller, Victor Ostromoukhov
 "Cascaded Sobol' Sampling", ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia), 40(6), pp. 274:1–274:13, December 2021
 December 2021
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 ``` 
