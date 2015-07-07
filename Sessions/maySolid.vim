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
badd +8 ConsState.cpp
badd +1 ConsState.h
badd +19 PrimState.cpp
badd +25 PrimState.h
badd +0 SolidSystem.cpp
badd +6 SolidSystem.h
badd +7 Makefile
badd +1 Solve1D.cpp
badd +0 Solve1D.h
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.cpp
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.h
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SquareMatrix.cpp
badd +1405 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Elastic2D.cpp
badd +1 Utils.cpp
badd +1 Utils.h
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticPrimState.cpp
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticPrimState.h
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Elastic1DUnigrid.cpp
badd +190 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SymmetricMatrix.cpp
badd +1 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticState.cpp
badd +45 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticState.h
badd +38 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/Romenski.cpp
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/Invariants.h
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SymmetricMatrix.h
badd +538 Elastic1D.cpp
badd +1 ElasticEOS.cpp
badd +21 ElasticEOS.h
badd +1 librarytest
badd +0 LibraryTest.cpp
badd +54 ElasticPrimState.h
badd +272 ~/Dropbox/2013-2014/Code/Solid/myCode/SimpleElastic/Solve1D.cpp
badd +4 fileRead.cpp
badd +14 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/ElasticEOS.h
badd +32 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticEOS/Romenski.h
badd +53 ElasticState.cpp
badd +43 ElasticState.h
badd +498 ~/Libraries/eigen/Eigen/src/Core/PermutationMatrix.h
badd +184 /usr/include/c++/4.6/bits/stl_bvector.h
badd +9 ./../SimpleElastic/ConsState.cpp
badd +3 ElasticPrimState1.h
badd +120 ElasticPrimState.cpp
args ConsState.cpp
edit Makefile
set splitbelow splitright
wincmd _ | wincmd |
split
1wincmd k
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe '1resize ' . ((&lines * 27 + 29) / 58)
exe '2resize ' . ((&lines * 27 + 29) / 58)
argglobal
let s:l = 7 - ((6 * winheight(0) + 13) / 27)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
7
normal! 0
wincmd w
argglobal
enew
wincmd w
exe '1resize ' . ((&lines * 27 + 29) / 58)
exe '2resize ' . ((&lines * 27 + 29) / 58)
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
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 6 - ((5 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
6
normal! 033l
wincmd w
argglobal
edit SolidSystem.cpp
let s:l = 44 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
44
normal! 053l
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
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
exe '1resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe '2resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 49 - ((21 * winheight(0) + 27) / 54)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
49
normal! 0
wincmd w
argglobal
edit ElasticState.cpp
let s:l = 115 - ((31 * winheight(0) + 27) / 54)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
115
normal! 09l
wincmd w
exe '1resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe '2resize ' . ((&lines * 54 + 29) / 58)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
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
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 57 - ((41 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
57
normal! 019l
wincmd w
argglobal
edit ElasticPrimState.cpp
let s:l = 49 - ((34 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
49
normal! 072l
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
tabedit ElasticState.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 1 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
edit ElasticState.h
let s:l = 43 - ((33 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
43
normal! 063l
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
tabedit Elastic1D.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 1 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
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
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 20 - ((19 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
20
normal! 033l
wincmd w
argglobal
edit ElasticEOS.cpp
let s:l = 1 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
tabedit LibraryTest.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 51 - ((33 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
51
normal! 014l
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Elastic1DUnigrid.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 1 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 105 - ((27 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
105
normal! 02l
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/System.cpp
let s:l = 497 - ((2 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
497
normal! 02l
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SquareMatrix.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 2 - ((1 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SymmetricMatrix.h
let s:l = 1 - ((0 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
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
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 296 - ((47 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
296
normal! 02l
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticPrimState.h
let s:l = 91 - ((38 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
91
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticState.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
argglobal
let s:l = 103 - ((51 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
103
normal! 02l
wincmd w
argglobal
edit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/ElasticState.h
let s:l = 9 - ((8 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
9
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 201)
exe 'vert 2resize ' . ((&columns * 100 + 100) / 201)
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
