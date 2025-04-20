#!/bin/bash

export MRB_BOT=$(pwd)
cd $MRB_TOP

mrbsetenv
mrb i -j10 

cd $MRB_BOT