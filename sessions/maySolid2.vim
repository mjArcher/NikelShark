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
badd +3 SolidSystem.cpp
badd +61 SolidSystem.h
badd +1 Makefile
badd +1 Solve1D.cpp
badd +1 Solve1D.h
badd +1088 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.cpp
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
badd +533 Elastic1D.cpp
badd +61 ElasticEOS.cpp
badd +63 ElasticEOS.h
badd +1 librarytest
badd +1 LibraryTest.cpp
badd +38 ElasticPrimState.h
badd +318 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Solve1D.cpp
badd +4 fileRead.cpp
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/ElasticEOS.h
badd +32 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/Romenski.h
badd +73 ElasticState.cpp
badd +49 ElasticState.h
badd +498 ~/Libraries/eigen/Eigen/src/Core/PermutationMatrix.h
badd +184 /usr/include/c++/4.6/bits/stl_bvector.h
badd +237 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/ConsState.cpp
badd +3 ElasticPrimState1.h
badd +103 ElasticPrimState.cpp
badd +32 /usr/include/c++/4.6/bits/c++0x_warning.h
badd +40 ~/Libraries/eigen-dev/unsupported/Eigen/CXX11/src/Core/util/CXX11Workarounds.h
badd +1 SolidUnitTests.cpp
badd +4 SolidEntroptBased.cpp
badd +89 /local/data/public/ma595/cns_amr/CNS_AMR_Multimaterial/modules/Problem/SolidEntropyBased/SolidEntropyBasedTest.C
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor3.h
badd +31 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor3.cpp
badd +1 Elastic1D.o:/lsc/zeushome/ma595/Dropbox/2013-2014/Code/Solid/myCode/dev/Elastic1D.cpp
badd +1 MA
badd +7 test.cpp
badd +1 ElasticState2.h
badd +4 ./locateFlags.sh
badd +0 ./Scripts/Solid1DAll.sh
badd +33 Makefile_1
badd +40 README
badd +1 ../../KevinElastic/Utils/Tensor4.cpp
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor4.h
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
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 31 + 104) / 209)
exe '2resize ' . ((&lines * 48 + 29) / 58)
exe 'vert 2resize ' . ((&columns * 88 + 104) / 209)
exe '3resize ' . ((&lines * 48 + 29) / 58)
exe 'vert 3resize ' . ((&columns * 88 + 104) / 209)
exe '4resize ' . ((&lines * 6 + 29) / 58)
exe 'vert 4resize ' . ((&columns * 177 + 104) / 209)
argglobal
enew
file NERD_tree_1
wincmd w
argglobal
let s:l = 58 - ((35 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
58
normal! 0
wincmd w
argglobal
edit Makefile
let s:l = 32 - ((29 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
32
normal! 0
wincmd w
argglobal
enew
wincmd w
4wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 104) / 209)
exe '2resize ' . ((&lines * 48 + 29) / 58)
exe 'vert 2resize ' . ((&columns * 88 + 104) / 209)
exe '3resize ' . ((&lines * 48 + 29) / 58)
exe 'vert 3resize ' . ((&columns * 88 + 104) / 209)
exe '4resize ' . ((&lines * 6 + 29) / 58)
exe 'vert 4resize ' . ((&columns * 177 + 104) / 209)
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
exe '1resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 1resize ' . ((&columns * 31 + 104) / 209)
exe '2resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 2resize ' . ((&columns * 177 + 104) / 209)
argglobal
enew
file NERD_tree_2
wincmd w
argglobal
let s:l = 731 - ((48 * winheight(0) + 27) / 54)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
731
normal! 0
wincmd w
4wincmd w
exe '1resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 1resize ' . ((&columns * 31 + 104) / 209)
exe '2resize ' . ((&lines * 54 + 29) / 58)
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
let s:l = 88 - ((50 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
88
normal! 0
wincmd w
argglobal
edit SolidSystem.cpp
let s:l = 133 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
133
normal! 0
wincmd w
4wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
tabedit ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/ConsState.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe '1resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe '2resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
argglobal
let s:l = 133 - ((0 * winheight(0) + 27) / 54)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
133
normal! 0
wincmd w
argglobal
edit ElasticEOS.cpp
let s:l = 80 - ((43 * winheight(0) + 27) / 54)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
80
normal! 04l
wincmd w
4wincmd w
exe '1resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe '2resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor4.h
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
let s:l = 4 - ((3 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
4
normal! 0
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Tensor3.h
let s:l = 1 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
4wincmd w
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
let s:l = 22 - ((21 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
22
normal! 02l
wincmd w
argglobal
edit ElasticState.cpp
let s:l = 97 - ((42 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
97
normal! 06l
wincmd w
4wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 108 + 104) / 209)
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
let s:l = 55 - ((39 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
55
normal! 01l
wincmd w
argglobal
edit ElasticPrimState.cpp
let s:l = 46 - ((31 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
46
normal! 0
wincmd w
4wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
tabedit LibraryTest.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 82 - ((30 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
82
normal! 0
4wincmd w
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Elastic1DUnigrid.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 154 - ((33 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
154
normal! 0
4wincmd w
tabnext 1
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
