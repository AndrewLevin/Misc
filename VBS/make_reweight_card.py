print "#******************************************************************"
print "#                       Reweight Module                           *"
print "#******************************************************************"
print "# launch"
print "#* Use the set command to specify the new set of parameter"
print "#* Or specify a path to a valid param_card/banner"
print "#* Example of valid command:"
print "#*     set aewm1 137"
print "#*     ~/param_card.dat"
print "#*"
print "#* Note:"
print "#*   1) the value of alphas will be used from the event"
print "#*      so the value of the param_card is not taken into account."
print "#*   2) It is dangerous to change a mass of any particle."
print ""
print ""
print "#* If you want to compute the weight for more than one hyppothesis"
print "#* you need first to uncomment the following line:"
print "# launch"
print "# and then use the set command to specify your parameter."
print "# All modification will start from the ORIGINAL card not from the"
print "# last define one."
print "#* You can have as many weight as you want."

for i in range(0,9):
    for j in range(0,9):
        if i == 0 and j == 0:
            continue
        print ""
        print "launch"
        print "        set anoinputs 1 " +str(i*5/10.) + "00000e-11"
        print "        set anoinputs 2 " + str(j*5/10.)+ "00000e-11"

exit(0)
