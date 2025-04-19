cd $MRB_TOP

mrbsetenv
mrb i -j10 #VERBOSE=1 #&> build_log.txt

#mv build_log.txt /exp/uboone/app/users/nlane/production/KaonShortProduction04/srcs/ubana/ubana/searchingforstrangeness

cd $MRB_TOP/srcs/ubana/ubana/searchingforstrangeness