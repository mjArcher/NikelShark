let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Dropbox/2013-2014/Code/Solid/myCode/dev
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +1 ConsState.cpp
badd +1 ConsState.h
badd +19 PrimState.cpp
badd +25 PrimState.h
badd +169 SolidSystem.cpp
badd +65 SolidSystem.h
badd +11 Makefile
badd +1 Solve1D.cpp
badd +1 Solve1D.h
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.cpp
badd +105 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.h
badd +2 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SquareMatrix.cpp
badd +1405 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Elastic2D.cpp
badd +23 Utils.cpp
badd +14 Utils.h
badd +231 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticPrimState.cpp
badd +13 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticPrimState.h
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Elastic1DUnigrid.cpp
badd +190 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SymmetricMatrix.cpp
badd +43 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticState.cpp
badd +28 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticState.h
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/Romenski.cpp
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Invariants.h
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SymmetricMatrix.h
badd +1 Elastic1D.cpp
badd +61 ElasticEOS.cpp
badd +63 ElasticEOS.h
badd +1 librarytest
badd +1 LibraryTest.cpp
badd +1 ElasticPrimState.h
badd +318 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Solve1D.cpp
badd +4 fileRead.cpp
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/ElasticEOS.h
badd +32 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/Romenski.h
badd +73 ElasticState.cpp
badd +1 ElasticState.h
badd +498 ~/Libraries/eigen/Eigen/src/Core/PermutationMatrix.h
badd +184 /usr/include/c++/4.6/bits/stl_bvector.h
badd +1 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/ConsState.cpp
badd +3 ElasticPrimState1.h
badd +103 ElasticPrimState.cpp
badd +32 /usr/include/c++/4.6/bits/c++0x_warning.h
badd +40 ~/Libraries/eigen-dev/unsupported/Eigen/CXX11/src/Core/util/CXX11Workarounds.h
badd +1 SolidUnitTests.cpp
badd +4 SolidEntroptBased.cpp
badd +89 /local/data/public/ma595/cns_amr/CNS_AMR_Multimaterial/modules/Problem/SolidEntropyBased/SolidEntropyBasedTest.C
badd +106 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor3.h
badd +31 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor3.cpp
badd +1 Elastic1D.o:/lsc/zeushome/ma595/Dropbox/2013-2014/Code/Solid/myCode/dev/Elastic1D.cpp
badd +1 MA
badd +7 test.cpp
badd +1 ElasticState2.h
badd +4 ./locateFlags.sh
badd +1 ./Scripts/Solid1DAll.sh
badd +46 Makefile_1
badd +40 README
badd +1 ../../KevinElastic/Utils/Tensor4.cpp
badd +4 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor4.h
badd +1 Tensor3.h
badd +6 SymmetricMatrix.h
badd +1 SymmetricMatrix.cpp
badd +4 SquareMatrix.h
badd +6 Tensor4.h
badd +1 Tensor.h
badd +65 ~/.vimrc
args ConsState.cpp
edit SolidSystem.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 46 + 30) / 60)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 10 + 30) / 60)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 214 - ((33 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
214
normal! 011l
wincmd w
argglobal
edit SolidSystem.h
let s:l = 1 - ((0 * winheight(0) + 23) / 46)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
enew
wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 46 + 30) / 60)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 10 + 30) / 60)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
tabedit Elastic1D.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 31 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 177 + 104) / 209)
argglobal
enew
file NERD_tree_2
wincmd w
argglobal
let s:l = 733 - ((53 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
733
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 177 + 104) / 209)
tabedit SolidSystem.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
argglobal
let s:l = 65 - ((9 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
65
normal! 0
wincmd w
argglobal
edit SolidSystem.cpp
let s:l = 215 - ((37 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
215
normal! 056l
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
tabedit ElasticEOS.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe '1resize ' . ((&lines * 56 + 30) / 60)
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 56 + 30) / 60)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 63 - ((54 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
63
normal! 0
wincmd w
argglobal
edit ElasticEOS.cpp
let s:l = 138 - ((37 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
138
normal! 0
wincmd w
exe '1resize ' . ((&lines * 56 + 30) / 60)
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 56 + 30) / 60)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
tabedit ElasticPrimState.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 33 - ((32 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
33
normal! 082l
wincmd w
argglobal
edit ElasticPrimState.cpp
let s:l = 140 - ((31 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
140
normal! 042l
wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
tabedit ElasticState.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
argglobal
let s:l = 22 - ((21 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
22
normal! 02l
wincmd w
argglobal
edit ElasticState.cpp
let s:l = 154 - ((52 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
154
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
tabedit LibraryTest.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 79 - ((52 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
79
normal! 026l
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Elastic1DUnigrid.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 154 - ((34 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
154
normal! 0
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 863 - ((28 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
863
normal! 024l
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticPrimState.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 273 - ((50 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
273
normal! 049l
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticPrimState.cpp
let s:l = 259 - ((50 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
259
normal! 024l
wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/Romenski.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 158 - ((45 * winheight(0) + 28) / 57)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
158
normal! 02l
tabnext 11
if exists('s:wipebuf')
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToOI
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
