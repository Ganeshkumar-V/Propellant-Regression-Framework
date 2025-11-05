#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | www.openfoam.com
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2011-2016 OpenFOAM Foundation
#     Copyright (C) 2016-2021 OpenCFD Ltd.
#     Copyright (C) 2023 Ganeshkumar V and Dilip S Sundaram.
#------------------------------------------------------------------------------
# License
#     This file is not a part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# File
#     etc/bashrc
#
# Description
#     The OpenFOAM environment for POSIX shell (eg, bash,dash,zsh,...).
#     Source manually or from the ~/.bashrc or ~/.profile files.
#
#------------------------------------------------------------------------------
# (advanced / legacy)
export PRF_PROJECT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export FOAM_PRF_APPBIN="$PRF_PROJECT_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_PRF_LIBBIN="$PRF_PROJECT_DIR/platforms/$WM_OPTIONS/lib"
export PATH=$FOAM_PRF_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_PRF_LIBBIN:$LD_LIBRARY_PATH
chmod +x Allwmake Allwclean
#------------------------------------------------------------------------------
