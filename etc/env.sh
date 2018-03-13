# for nwpy:
export PYTHONPATH=$PYTHONPATH:/home/Javier.Delgado/libs/nwpy/lib
# for Equation:
echo $PYTHONPATH | grep ":/home/Javier.Delgado/local/lib/python2.7/site-packages:" &> /dev/null
[[ $? != 0 ]] && export PYTHONPATH=$PYTHONPATH:/home/Javier.Delgado/local/lib/python2.7/site-packages
# pycane
export PYTHONPATH=/home/Javier.Delgado/apps/pycane_dist/trunk:$PYTHONPATH
