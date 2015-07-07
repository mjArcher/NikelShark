let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Code/Solid/myCode/dev
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +227 ElasticPrimState.cpp
badd +1 ./Elastic1D.cpp
badd +132 ElasticEOS.cpp
badd +214 SolidSystem.cpp
badd +1 ./../../KevinElastic/System.cpp
badd +1 ./../../KevinElastic/ElasticPrimState.cpp
badd +41 ElasticPrimState.h
badd +1 ./LibraryTest.cpp
badd +1 ./librarytest
badd +1 ./LibraryTest.sh
badd +13 ./Scripts/Solid1DAll.sh
badd +4 SquareTensor3.h
badd +10 ./TensorLibraryOld/Tensor3.cpp
badd +0 ./TensorLibraryOld/Tensor3.h
badd +0 SquareTensor3.cpp
badd +58 SolidSystem.h
badd +0 Makefile
badd +49 ~/Libraries/eigen/Eigen/src/Core/../plugins/CommonCwiseUnaryOps.h
badd +0 ./../../KevinElastic/Elastic1DUnigrid.cpp
badd +0 ElasticState.
badd +0 ElasticState.cpp
badd +33 ElasticEOS.h
badd +410 ~/Libraries/eigen/Eigen/src/Core/PlainObjectBase.h
args ElasticPrimState.cpp
edit SolidSystem.h
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
exe '2resize ' . ((&lines * 37 + 29) / 59)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 18 + 29) / 59)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 21 - ((0 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
21
normal! 066l
wincmd w
argglobal
edit ElasticEOS.cpp
let s:l = 163 - ((27 * winheight(0) + 18) / 37)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
163
normal! 0
wincmd w
argglobal
enew
wincmd w
2wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 37 + 29) / 59)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 18 + 29) / 59)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
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
let s:l = 41 - ((26 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
41
normal! 071l
wincmd w
argglobal
edit ElasticPrimState.cpp
let s:l = 228 - ((37 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
228
normal! 0
wincmd w
2wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
tabedit ./Elastic1D.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 1 - ((0 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
2wincmd w
tabedit ElasticEOS.h
set splitbelow splitright
wincmd _ | wincmd |
split
1wincmd k
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe '1resize ' . ((&lines * 0 + 29) / 59)
exe '2resize ' . ((&lines * 55 + 29) / 59)
argglobal
let s:l = 33 - ((21 * winheight(0) + 0) / 0)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
33
normal! 042l
wincmd w
argglobal
edit ElasticEOS.cpp
let s:l = 157 - ((48 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
157
normal! 01l
wincmd w
2wincmd w
exe '1resize ' . ((&lines * 0 + 29) / 59)
exe '2resize ' . ((&lines * 55 + 29) / 59)
tabedit SolidSystem.h
set splitbelow splitright
wincmd _ | wincmd |
split
1wincmd k
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe '1resize ' . ((&lines * 0 + 29) / 59)
exe '2resize ' . ((&lines * 55 + 29) / 59)
argglobal
let s:l = 71 - ((53 * winheight(0) + 0) / 0)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
71
normal! 079l
wincmd w
argglobal
edit SolidSystem.cpp
let s:l = 200 - ((19 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
200
normal! 0
wincmd w
2wincmd w
exe '1resize ' . ((&lines * 0 + 29) / 59)
exe '2resize ' . ((&lines * 55 + 29) / 59)
tabedit Makefile
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 39 - ((38 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
39
normal! 068l
2wincmd w
tabedit ./LibraryTest.sh
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
exe '2resize ' . ((&lines * 48 + 29) / 59)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 7 + 29) / 59)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 4 - ((3 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
4
normal! 0
wincmd w
argglobal
edit ./LibraryTest.cpp
let s:l = 94 - ((45 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
94
normal! 024l
wincmd w
argglobal
enew
wincmd w
2wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 48 + 29) / 59)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 7 + 29) / 59)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
tabedit ElasticState.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 123 - ((38 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
123
normal! 013l
2wincmd w
tabedit ./TensorLibraryOld/Tensor3.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd _ | wincmd |
split
1wincmd k
wincmd w
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe '1resize ' . ((&lines * 0 + 29) / 59)
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 55 + 29) / 59)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 55 + 29) / 59)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
exe '4resize ' . ((&lines * 0 + 29) / 59)
exe 'vert 4resize ' . ((&columns * 104 + 104) / 209)
argglobal
let s:l = 6 - ((4 * winheight(0) + 0) / 0)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
6
normal! 0
wincmd w
argglobal
edit ./TensorLibraryOld/Tensor3.h
let s:l = 135 - ((50 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
135
normal! 0
wincmd w
argglobal
edit SquareTensor3.cpp
let s:l = 31 - ((30 * winheight(0) + 27) / 55)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
31
normal! 0
wincmd w
argglobal
edit SquareTensor3.h
let s:l = 65 - ((49 * winheight(0) + 0) / 0)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
65
normal! 059l
wincmd w
2wincmd w
exe '1resize ' . ((&lines * 0 + 29) / 59)
exe 'vert 1resize ' . ((&columns * 104 + 104) / 209)
exe '2resize ' . ((&lines * 55 + 29) / 59)
exe 'vert 2resize ' . ((&columns * 104 + 104) / 209)
exe '3resize ' . ((&lines * 55 + 29) / 59)
exe 'vert 3resize ' . ((&columns * 104 + 104) / 209)
exe '4resize ' . ((&lines * 0 + 29) / 59)
exe 'vert 4resize ' . ((&columns * 104 + 104) / 209)
tabedit ./../../KevinElastic/System.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 408 - ((19 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
408
normal! 0
2wincmd w
tabedit ./../../KevinElastic/ElasticPrimState.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 315 - ((55 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
315
normal! 0
2wincmd w
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
