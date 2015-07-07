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
badd +0 ConsState.cpp
badd +0 ./ConsState.h
badd +0 ./PrimState.cpp
badd +0 PrimState.h
badd +0 ./SolidSystem.cpp
badd +0 SolidSystem.h
badd +0 ./Makefile
badd +0 ./Solve1D.cpp
badd +0 Solve1D.h
badd +0 ./../../KevinElastic/System.cpp
badd +0 ./../../KevinElastic/System.h
badd +0 ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SquareMatrix.cpp
badd +0 ./../../KevinElastic/Elastic2D.cpp
badd +0 Utils.cpp
badd +4 Utils.h
badd +0 ./../../KevinElastic/ElasticPrimState.cpp
args ConsState.cpp
edit ./ConsState.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 100) / 200)
exe 'vert 2resize ' . ((&columns * 99 + 100) / 200)
argglobal
let s:l = 53 - ((52 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
53
normal! 0
wincmd w
argglobal
edit ConsState.cpp
let s:l = 1 - ((0 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 200)
exe 'vert 2resize ' . ((&columns * 99 + 100) / 200)
tabedit PrimState.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 100 + 100) / 200)
exe 'vert 2resize ' . ((&columns * 99 + 100) / 200)
argglobal
let s:l = 1 - ((0 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
edit ./PrimState.cpp
let s:l = 1 - ((0 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 200)
exe 'vert 2resize ' . ((&columns * 99 + 100) / 200)
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
exe 'vert 1resize ' . ((&columns * 100 + 100) / 200)
exe 'vert 2resize ' . ((&columns * 99 + 100) / 200)
argglobal
let s:l = 11 - ((10 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
11
normal! 0
wincmd w
argglobal
edit ./SolidSystem.cpp
let s:l = 3 - ((2 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
3
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 100 + 100) / 200)
exe 'vert 2resize ' . ((&columns * 99 + 100) / 200)
tabedit ./Solve1D.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 16 - ((15 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
16
normal! 0
tabedit ./Makefile
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 10 - ((9 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
10
normal! 069l
tabedit ./../../KevinElastic/System.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 763 - ((24 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
763
normal! 08l
tabedit ./../../KevinElastic/ElasticPrimState.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 11 - ((10 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
11
normal! 078l
tabedit ./../../KevinElastic/System.h
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 1 - ((0 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
tabedit ~/Dropbox/2013-2014/Code/Solid/KevinElastic/Utils/SquareMatrix.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 1 - ((0 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
tabedit ./../../KevinElastic/Elastic2D.cpp
set splitbelow splitright
set nosplitbelow
set nosplitright
wincmd t
set winheight=1 winwidth=1
argglobal
let s:l = 1 - ((0 * winheight(0) + 26) / 53)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
tabnext 6
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
