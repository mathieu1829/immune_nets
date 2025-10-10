let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd ~/Documents/python/immune_nets
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
let s:shortmess_save = &shortmess
if &shortmess =~ 'A'
  set shortmess=aoOA
else
  set shortmess=aoO
endif
badd +1 src/orm/models/repertoire.py
badd +1 src/orm/models/dataset.py
badd +0 tests/test_orm.py
badd +23 src/orm/models/clonotypeData.py
argglobal
%argdel
edit src/orm/models/dataset.py
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
2wincmd h
wincmd _ | wincmd |
split
1wincmd k
wincmd w
wincmd w
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe '1resize ' . ((&lines * 37 + 38) / 76)
exe 'vert 1resize ' . ((&columns * 126 + 190) / 381)
exe '2resize ' . ((&lines * 36 + 38) / 76)
exe 'vert 2resize ' . ((&columns * 126 + 190) / 381)
exe 'vert 3resize ' . ((&columns * 126 + 190) / 381)
exe '4resize ' . ((&lines * 37 + 38) / 76)
exe 'vert 4resize ' . ((&columns * 127 + 190) / 381)
exe '5resize ' . ((&lines * 36 + 38) / 76)
exe 'vert 5resize ' . ((&columns * 127 + 190) / 381)
argglobal
balt src/orm/models/repertoire.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 18) / 36)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
argglobal
if bufexists(fnamemodify("tests/test_orm.py", ":p")) | buffer tests/test_orm.py | else | edit tests/test_orm.py | endif
if &buftype ==# 'terminal'
  silent file tests/test_orm.py
endif
balt src/orm/models/dataset.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
argglobal
if bufexists(fnamemodify("src/orm/models/repertoire.py", ":p")) | buffer src/orm/models/repertoire.py | else | edit src/orm/models/repertoire.py | endif
if &buftype ==# 'terminal'
  silent file src/orm/models/repertoire.py
endif
balt src/orm/models/repertoire.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 36) / 73)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
argglobal
if bufexists(fnamemodify("src/orm/models/repertoire.py", ":p")) | buffer src/orm/models/repertoire.py | else | edit src/orm/models/repertoire.py | endif
if &buftype ==# 'terminal'
  silent file src/orm/models/repertoire.py
endif
balt src/orm/models/repertoire.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 1 - ((0 * winheight(0) + 18) / 36)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 1
normal! 0
wincmd w
argglobal
if bufexists(fnamemodify("src/orm/models/clonotypeData.py", ":p")) | buffer src/orm/models/clonotypeData.py | else | edit src/orm/models/clonotypeData.py | endif
if &buftype ==# 'terminal'
  silent file src/orm/models/clonotypeData.py
endif
balt src/orm/models/repertoire.py
setlocal foldmethod=manual
setlocal foldexpr=0
setlocal foldmarker={{{,}}}
setlocal foldignore=#
setlocal foldlevel=0
setlocal foldminlines=1
setlocal foldnestmax=20
setlocal foldenable
silent! normal! zE
let &fdl = &fdl
let s:l = 23 - ((22 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 23
normal! 050|
wincmd w
5wincmd w
exe '1resize ' . ((&lines * 37 + 38) / 76)
exe 'vert 1resize ' . ((&columns * 126 + 190) / 381)
exe '2resize ' . ((&lines * 36 + 38) / 76)
exe 'vert 2resize ' . ((&columns * 126 + 190) / 381)
exe 'vert 3resize ' . ((&columns * 126 + 190) / 381)
exe '4resize ' . ((&lines * 37 + 38) / 76)
exe 'vert 4resize ' . ((&columns * 127 + 190) / 381)
exe '5resize ' . ((&lines * 36 + 38) / 76)
exe 'vert 5resize ' . ((&columns * 127 + 190) / 381)
tabnext 1
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20
let &shortmess = s:shortmess_save
let &winminheight = s:save_winminheight
let &winminwidth = s:save_winminwidth
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
set hlsearch
nohlsearch
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
