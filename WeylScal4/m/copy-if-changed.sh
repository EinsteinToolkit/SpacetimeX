#! /bin/bash

set -e

src=$1
dst=$2

# Copy tree $src to tree $dst


if test -z "$src" || test -z "$dst" || test "$src" = "$dst"; then
    echo "Usage: $0 <src> <dst>"
    exit 1
fi



if test -f $src; then
    
    # Both $src and $dst must be files.
    
    test -f $dst || ! test -e $dst || exit 3
    
    # Copy file if it differs
    srcfile=$src
    dstfile=$dst
    if cmp -s $srcfile $dstfile; then
        : # unchanged; do nothing
    else
        echo cp $srcfile
        cp $srcfile $dstfile
    fi
    
fi



if test -d $src; then
    
    # Both $src and $dst must be directories.  $dst is created if it
    # does not exist
    
    # All files in the source tree are checked; if they already exist
    # in the destination tree and are identical, they are ignored,
    # otherwise they are copied.  Missing directories are created.
    
    # All files in the destination tree are checked; if they do not
    # exist in the source tree, they are deleted.
    
    test -e $dst || mkdir -p $dst
    test -d $dst || exit 3

    # find inode number to identify the src directory later on
    srcinode=$(ls -id $src | awk '{print $1}')

    # Create all directories
    for dir in $(cd $src && find . -type d); do
        dstdir=$dst/$dir
        if test -d $dstdir; then
            : # directory exists; do nothing
        else
            echo mkdir $dstdir
            mkdir -p $dstdir
        fi
    done
    
    # Delete directories which do not exist
    for dir in $(cd $dst && find . \( -name 'CVS' -o -name '.svn' -o -name '.git' -o -name '.hg' -o -name '_darcs' -o -inum $srcinode \) -prune -o -type d -print); do
        srcdir=$src/$dir
        dstdir=$dst/$dir
        if test -d $srcdir; then
            : # directory exists; do nothing
        else
            echo rm -rf $dstdir
            rm -rf $dstdir
        fi
    done
    
    # Copy files that differ
    for file in $(cd $src && find . -type f); do
        srcfile=$src/$file
        dstfile=$dst/$file
        if cmp -s $srcfile $dstfile; then
            : # unchanged; do nothing
        else
            echo cp $srcfile
            cp $srcfile $dstfile
        fi
    done
    
    # Delete files which do not exist
    for file in $(cd $dst && find . \( -name 'CVS' -o -name '.svn' -o -name '.git' -o -name '.hg' -o -name '_darcs' -o -inum $srcinode \) -prune -o -type f -print); do
        srcfile=$src/$file
        dstfile=$dst/$file
        if test -e $srcfile; then
            : # file exists; do nothing
        else
            echo rm $dstfile
            rm $dstfile
        fi
    done
    
fi
