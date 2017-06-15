__author__ = 'chad'

import re, sys

#This script summarizes the results of an Kick-R search.
#A usage example is 'python kicksum.py 250 > results', which will create a file 'results' that summarizes your search
#based on 250 candidate structures.


def main():

    jobs_searched = int(sys.argv[1])
    jobs_run = 0
    jobs_converged = 0
    converged_list = []

    print "Summary of Output Files Checked"
    print"************************************************************\n"

    for i in range(1, jobs_searched + 1):
        try:
            read_in = open(str(i) + '.kick.out', 'r')
            output = ''.join([lines for lines in read_in])
            output = output.replace('\n', '').replace(' ', '')

            if re.search('HF=(-?\d+.\d+)', output) is not None:
                termination_status = ' Converged: E = '
                energy = re.search('HF=(-?\d+.\d+)', output)
                print str(i) + termination_status + energy.group(1)
                jobs_run = jobs_run + 1
                jobs_converged = jobs_converged + 1
                converged_list.append([i, float(energy.group(1))])
            else:
                termination_status = ' Failed'
                print str(i) + termination_status
                jobs_run = jobs_run + 1

        except IOError:
            print str(i) + ' No Output'

        i = i + 1

    convergence_percent = (jobs_converged/float(jobs_run))*100
    converged_list.sort(key=lambda x: x[1])
    unique_structures = []
    unique_structures.append(converged_list[0])

    for i in range(1,len(converged_list)):
        if abs(unique_structures[len(unique_structures)-1][1] - converged_list[i][1]) > 0.00001:
            unique_structures.append(converged_list[i])

    lowest_energy = unique_structures[0][1]

    print"\n************************************************************"
    print "\n" + "Jobs Searched: " + str(jobs_searched)
    print "Jobs Run: " + str(jobs_run)
    print "Jobs Converged: " + str(jobs_converged)
    print "Convergence Percentage: %" + str(convergence_percent)
    print "Unique Structures: " + str(len(unique_structures)) + "\n"
    print"************************************************************\n"
    print "Unique Structures\n"

    for i in range(len(unique_structures)):
        relative_energy = unique_structures[i][1] - lowest_energy
        print "Job Number: " + str(unique_structures[i][0])
        print "Energy: " + str(unique_structures[i][1]) + "   (Hartrees)"
        print "Relative Energy: +" + str(relative_energy*627.5095) + "   (kcal/mol)"
        print "Geometry: "
        read_in = open(str(unique_structures[i][0]) + '.kick.out', 'r')
        output = ''.join([lines for lines in read_in])
        output = output.replace('\n', '').replace(' ', '')
        templist = output.split("\\\\")
        templist = templist[3].split("\\")
        for i in range(len(templist)):
            templist[i] = templist[i].replace(",", "   ")
        for i in range(1, len(templist)):
            print templist[i]
        print "\n"




main()


