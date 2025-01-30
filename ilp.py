# ASSIGNMENT 3

import numpy as np
import math
def gomory(filename):
    """
    2 PHASE SIMPLEX METHOD IMPLEMENTATION:

    PHASE I:
    1. Convert the problem to standard form
    2. By multiplying by (-1), make sure that b>=0
    3. Introduce the artificial variables y1 upto ym and apply the simplex method to the auxiliary problem
    4. If optimal cost of auxiliary problem is positive, problem is infeasible
    5. If optimal cost is zero, a feasible soln to the original problem has been found
    6. If no artificial variables in the final basis, articificial variables and corr. columns eliminated
    7. If lth basic variable is an artificial one, examine lth entry of columns B-1Aj. If all of these are zero, lth row represents a redundant constraint
    8. Otherwise, apply a change of basis 
    9. Repeat this operation to drive out all artificial variables

    PHASE II:
    1. let the final basis and tableau from Phase I be the initial basis and tableau for Phase II
    2. Compute the reduced costs of all variables for this initial basis, using the cost coefficients of the original problem
    3. Apply the simplex method to the original problem 
    """

    # BASED ON THE ABOVE PROCEDURE, WE DEFINE THE FOLLOWING HELPER FUNCTIONS:
    """
    -- 1. convert_std()    # converts to standard form
    -- 2. be_positive()    # (pun) but also makes sure that b is positive 
    -- 3. var_art()        # introduces artificial variables while also checking if the slack variables introduced earlier suffice
    -- 4. simplex_method() # full tableau implementation of the simplex method with the initial bfs st all original variables are zero and artificial variables are equal to b
    -- 5. is_feasible()    # checks the value of optimal cost to see if the problem is feasible 
    -- 6. basic_art()      # checks to see if there are any artificial variables in the basis and returns the index of the first art basic var it encounters
    -- 7. check_zero()     # checks if the artificial variable corresponds to a redundant constraint 
    -- 8. change_basis()   # given the pivot row and column, it performs required operations for entry and exit of variables from the basis  
    10. re_reduce()     # (reuse, reduce, re_reduce) recomputes the reduced costs for the problem in phase II
    """
    # We also have our main function called phase_2_simplex() which combines all of these to execute 2-phase simplex 

    # the following function reads the input file 
    def read_file(filename):
        A=[]
        f= open(filename,'r+')
        l1=list(map(int,f.readline().split()))
        b=list(map(int,f.readline().split()))
        c=list(map(int,f.readline().split()))
        for i in range (0,len(b)):
            A.append(list(map(int,f.readline().split())))
        f.close()
        c=np.array(c)
        c=c*(-1)
        out=(l1[0],l1[1],c,A,b)
        return out
    #print(read_file("input.txt"))


    def num_art(b):
        count = 0
        # this function counts the number of artificial variables that will be needed
        for i in range(0, len(b)):
            if b[i]<0:
                count = count +1
        return count

    def be_positive(table, art):
        # here, art represents the list of rows which have to be multiplied by -1 
        for j in range(0, len(art)):
            i = art[j]
            table[i,:] = (-1)*table[i,:]
        return table 

    def var_art(table, art_add, n, m, basis):
        for j in range (0, len(art_add)):
            i = art_add[j]
            table[i, n+m+j+1] = 1
            basis.append(n+m+j+1)
        return (table, basis) 

    def convert_std(A, b, c, n, m):
        # we are given a matrix A and 2 vectors b and c along with n variables and m constraints 
        # the function returns a tableau as np array 
        a = num_art(b)
        #print(a)
        rows = m+1
        columns = n+m+a+1 # an additional slack variable is added for every constraint here 
        table = np.zeros((rows, columns))
        for i in range(1,rows):
            table[i,0] = b[i-1]
        #table[:,0] = b     # we set the first column as b
        for i in range(1,m+1):
            for j in range (1, n+1):
                table[i,j] = A[i-1][j-1]    # we add the coeff matrix to the tableau as well 
        # we will now have to check the artificial and slack variables and see whether b values are -ve
        basis = []   # this will be a list of basis variables 
        art_add = []
        for i in range(0, m):
            if b[i] >=0:
                table[i+1,n+i+1] = 1
                basis.append(n+i+1) # adding the index of the column to the basis vectors list
            else:
                table[i+1,n+i+1] = 1
                # we will also have to multiply the entire row by -1, we also add both a slack variable and an artificial variable 
                art_add.append(i+1)   # adding to the list of indices corresponding to which artificial variables have to be added 
                #basis.append(n+i+1)
        # we now make all the b's positive 
        table = be_positive(table, art_add)
        #print(f'table after accounting for negative b: {table}')    # accounting hasn't been done 
        # we now add the artificial variables to the table 
        (table, basis) = var_art(table, art_add, n, m, basis)
        # after all of this has been done, we still need to take care of the zeroth row 
        # we know that the basis matrix initially is the identity matrix
        # the reduced costs are zero for all the basic variables, and for all the non basic variables, they are -(cB)T(A), where cB  = 1 because of the new function introduced 
        # so, we basically sum up the entries in the corresponding column 
        #print(table, basis)
        #print(f'c:{c}')
        #print
        for i in range(0, columns):
            if i not in basis:
                # if i is not in basis, then we add up the columns 
                sum = 0
                for j in range(1,rows):
                    sum = sum + table[j,i]
                #print(i)
                if i<=len(c):
                    table[0,i] = -sum
                else:
                    table[0,i] = -sum
        # hence, we have computed the reduced costs
        # we now set the total cost, which will be negative of the sum of b values 
        #table[0,0] = 0
        """
        cost = 0
        for i in range(1, rows):
            cost = cost + b[i-1]
        table[0,0] = -cost
        """
        # with this, our table is ready for simplex!
        #print(f'standardized table: {table}')
        art = []
        for i in range (0, len(art_add)):
            art.append(art_add[i]+n+m)
        return (table, basis, art)

    # TESTED AND WORKS!

    def is_feasible(table):
        if table[0,0] == 0:
            return True
        else:
            return False 
        
    def change_basis(table, piv_r, piv_c):
        # given the indices of the pivot row and column, it performs change of basis operations
        # variable corresponding to the row leaves the basis and that corresponding to the column enters it
        # we multiply the piv_row such that all elements in the pivot column except the pivot element become 0
        # the pivot element becomes 1
        pivot = table[piv_r,piv_c]  # this is the pivot element
        #print(f'pivot: {pivot}')
        piv_row = table[piv_r, :]
        piv_column = table[:, piv_c]
        piv_row = piv_row/pivot
        table[piv_r,:] = piv_row
        for i in range(0,len(piv_column)):
            if i!=piv_r:
                # row = table[i,:]
                coeff = table[i,:][piv_c]
                table[i,:] = table[i,:] - coeff*piv_row
        return table 

    def is_finished(table):
        # checks if the algorithm should terminate
        r = len(table[0,:])
        zero_r = table[0,:]
        i = 1
        while i<r:
            if zero_r[i] < 0:
                return False
            else:
                i = i+1
        return True 

    def find_piv_row(table, piv_c):
        # given the table and the pivot column, it returns the index of the pivot row 
        l_col = len(table[:,0])
        piv_col = table[:,piv_c]
        zero_col = table[:,0]
        # print(piv_col)
        # print(zero_col)
        # we take the ratio of ith entry of the zeroth column and the ith entry of the pivot column
        # we need the minimum of these ratios
        #print(piv_col)
        i = 1
        found = False
        while i < l_col and not found:
            if piv_col[i] > 0:
                min = (zero_col[i]/piv_col[i])
                piv_r = i
                found = True
            else:
                i = i+1
        if not found:
            piv_r = 'problem unbounded'
        for i in range(i+1,l_col):
            if piv_col[i] > 0:
                ratio = (zero_col[i]/piv_col[i])
                # print(ratio)
                if min > ratio:
                    min = ratio
                    piv_r = i
        return piv_r

    def simplex_method(table, basis):
        # given the table and the basis, it computes the optimal solution using the full tableau method 
        row = len(table[0,:])
        col = len(table[:,0])
        while not is_finished(table):
            # this means while all the reduced costs are not non negative
            # we now choose the j for which the reduced cost is negative
            j = 1
            pc = 1
            found = False
            while j < row and not found:
                if table[0,j] < 0:
                    pc = j
                    found = True
                else:
                    j = j + 1
            # we now find the pivot row
            pr = find_piv_row(table, pc) 
            #print(f'pc:{pc}')
            # now, variable corresponding to pivot row exits the basis and that corresponding to pivot column enters the basis
            #print(basis)
            #print(f'pr:{pr}')
            basis[pr-1] = pc
            table = change_basis(table, pr, pc)
            #print(f'simplex table: {table}')
        # once the optimal solution has been found, we return the table 
        return (table, basis) 

    # NOW, WE TEST THE HELPER FUNCTIONS FOR SIMPLEX METHOD:
    """
    1. is_feasible 
    2. change_basis
    3. find_piv_row
    4. simplex_method
    """
    # convert_std working fine 


    # THIS IS FOLLOWED BY DRIVING OUT THE ARTIFICIAL VARIABLES THAT MIGHT BE PRESENT IN THE BASIS 

    def basic_art(basis, art):
        # checks to see if there are any artificial variables in the basis
        # here, art is the list of indices of artificial variables
        # it also returns the index of the first artificial variable it encounters in the basis 
        for i in basis:
            for j in art:
                if i == j:
                    return (True, i) 
        return (False, None) 

    # print(basic_art([1, 2, 3], [3, 6, 7]))
    # THIS WORKS!

    def check_zero(table, i, art):
        # checks to see if the constraint is redundant 
        row_l = len(table[0,:]) 
        #print(f'row_l:{row_l}')
        for j in range(1, row_l):
            # we check to see if it is an artificial variable 
            if (j-1) not in art:
                if table[i,j] != 0:
                    return False
        return True 

    def find_non_basic(table, basic, art):
        # finds the variable on the original problem that is not a part of the basis 
        for i in range(1, len(table[0,:])):
            if i not in basic:
                if i not in art:
                    return i

    def find_next_row(curr, l, rem_rows):
        while curr < l:
            if curr in rem_rows:
                curr = curr + 1
            else:
                return curr


    def remake_tableau(table, art, rem_rows):
        # keep all the variables except the artificial variables 
        #print(art)
        l = len(table[:,0])
        l_col_n = (len(table[0,:])-len(art))
        l_row_n = len(table[:,0]) - len(rem_rows)
        #print(l_col_n, l_row_n)
        # print(art, rem_rows)    # why is 'art' empty?
        # print(l_col_n, l_row_n) # l_col_n is incorrect - should be 4, printing 7
        new_t = np.zeros((l_row_n, l_col_n))
        # l_row_n = 4, l_col_n = 6
        # print(new_t)
        # once the new table has been created, we can start assigning values accordingly - we go ROW WISE 
        row_num = 1 # this is the index of the row in the previous table that has to be assigned 
        for i in range(1, l_row_n):
            #print(row_num)
            for j in range(0, l_col_n):
                new_t[i,j] = table[row_num, j]
            row_num = find_next_row(row_num+1, l, rem_rows)
            #print(f'new table: {new_t}')
            #print(f'j:{j}')
        return new_t

    def re_reduce(table, c, basis, n):
        # calculates reduced cost again
        cb = []
        row_l = len(table[0,:])
        col_l = len(table[:,0])
        #print(f'r,c:{row_l, col_l}')
        #print(f'n:{n}')
        for i in range(1,row_l):
            if i in basis:
                #print(i)
                if i <= len(c):
                    cb.append(c[i-1])
                else:
                    cb.append(0)
        #print(cb)
        # we now have cb
        cb_f = cb.copy()
        basis_n = basis.copy()
        # this cb_f is basically the cost vector values for basis vectors in the order in which they appear in the table for effective dot product calculation
        # basically find pos of element with minimum index and put the new element there 
        for i in range (0, len(basis)):
            minimum = min(basis_n)
            ind = basis.index(minimum)
            basis_n[ind] = 2*row_l
            cb_f[ind] = cb[i]
        #print(f'final basis: {cb_f}')
        for i in range(0, row_l):
            if i in basis:
                table[0,i] = 0
            else:
                if i == 0:
                    in_c = 0
                elif i <= n:
                    in_c = c[i-1]
                else:
                    in_c = 0
                sum = 0
                col = table[:,i]
                for j in range(1, col_l):
                    #print(f'cb, col:{cb_f[j-1], col[j]}') # got the problem here - not being multiplied by respective positions
                    #print(f'prod = {cb_f[j-1]*col[j]}')
                    sum = sum + cb_f[j-1]*col[j]
                    #print(f'sum = {sum, j}')
                #print(sum)
                cost = (in_c - sum)
                table[0,i] = cost
            #print(table)
        return table 

    def phase_2_simplex(A, b, c, n, m):
        #print(A,b,c,n,m)
        (table, basis, art) = convert_std(A,b,c,n,m)
        #print(f'initial basis: {basis}')
        #print(table)
        #print(table)
        # after converting to standard form and forming the tableau, we perform simplex method on it
        (table, basis) = simplex_method(table, basis)
        # now, we need to drive out the artificial variables 
        #print(f'table after auxiliary problem solved: {table}')
        rem_rows = []
        # print(f"art = {art}")   # convert_std is not returning artificial variables hi
        #print(basic_art(basis, art))
        while basic_art(basis, art)[0]:
            #print(f"art = {art}")   # this is incorrect - should be [9 10 11] but is [1 2 3]
            #print(f"basis = {basis}")
            # while there are artificial variables in the basis,
            #print(f'index: {basis.index(basic_art(basis, art)[1])}')
            piv_row = basis.index(basic_art(basis, art)[1])+1
            if not check_zero(table, piv_row, art):
                piv_col = find_non_basic(table, basis, art)
                #print(f'piv_col:{piv_col}')
                # we now use change of basis 
                table = change_basis(table, piv_row, piv_col)
                # we update the list accordingly 
                basis[piv_row-1] = piv_col
                #print(f'new_basis: {basis}')
            else:
                basis.remove(piv_row)   # since this is now a redundant constraint 
                rem_rows.append(piv_row)
        #print(f'table after the artificial variables have been driven out:{table}')
        # By the end of the while loop, the artificial variables would have been driven out of the basis 
        # we will now have to remake the tableau and re calculate the reduced costs
        # print("table huehue")
        # print(table)
        #print(f'removed rows: {rem_rows}')
        #print(table)
        table = remake_tableau(table, art, rem_rows)    # the issue lies here
        #print(f'table after the tableau has been remade: {table}')
        #print(table)
        # now we need to recalculate the reduced costs for this new table 
        # reduced cost for all basic variables is 0
        # reduced cost for non basic variables - (cj - (cB)T(B-1)(Aj))
        table = re_reduce(table, c, basis, n)
        #print(f'Table after the costs have been recalculated: {table}')
        (table, basis) = simplex_method(table, basis)
        #print(table) # this is the problem - re_reduce has not worked properly
        return(table, basis)

    def can_be_improved(tableau):
        for i in range(0,tableau.shape[0]):
            if tableau[i,0]<0 and i!=0:
                return [1,i]
        return [0,0]

    def get_pivot_position(tableau):
        res=can_be_improved(tableau)
        if res[0]==1:
            flag=0
            for i in range (1,tableau.shape[1]):
                if tableau[res[1]][i]<0:
                    min1=[tableau[0][i]/abs(tableau[res[1]][i]),i,res[1]]
                    if min1!=0:
                        if flag==0:
                            min=min1
                            flag=1
                        if min1<min:
                            min=min1
        if flag==1:
            return min
        else:
            return -1
        
    def modify_matrix(tableau,l):
        min=get_pivot_position(tableau)
        '''for i in range (len(l)):
            if l[i]==min[2]:
                l[i]=min[1]'''
        l[min[2]-1]=min[1]
        #print(l)
        for i in range (0,tableau.shape[0]):
            factor=tableau[i][min[1]]/tableau[min[2]][min[1]]
            for j in range (0,tableau.shape[1]):
                if i!=min[2]:
                    tableau[i][j]=tableau[i][j]-factor*tableau[min[2]][j]
        y=tableau[min[2]][min[1]]
        for j in range (0,tableau.shape[1]):
            x=tableau[min[2]][j]/y
            tableau[min[2]][j]=x
        return (tableau,l)

    def dual_simplex(tableau,l):
        res=can_be_improved(tableau)
        while res[0]!=0:
            (tableau,l)=modify_matrix(tableau,l)
            res=can_be_improved(tableau)
        return (tableau,l)

    #gomory cut algorithm
    def gomory_constraint(tableau):
        l=[]
        flag=0
        for i in range(0,tableau.shape[0]):
            if flag==0 and int(tableau[i][0])!=tableau[i][0]:
                for j in range(tableau.shape[1]):
                    fraction=tableau[i][j]%1
                    l.append(fraction)
                if l==[0]*tableau.shape[1]:
                    l=[]
                elif l!=[0]*tableau.shape[1]:
                    flag=1
            if flag==1:
                break
        return l

    def gomory_tableau(tableau,n):
        n.append(tableau.shape[1])
        l=gomory_constraint(tableau)
        l=np.array(l)
        l=l*-1
        tableau=np.row_stack((tableau,l))
        m=[0]*tableau.shape[0]
        m[-1]=1
        m=np.array(m)
        tableau=np.column_stack((tableau,m))
        return (tableau,n)

    '''def gomory_implementation(tableau):
        flag=0
        for i in range (tableau.shape[0]):
            if int(tableau[i][0])!=tableau[i][0] and flag==0:
                tableau=gomory_tableau(tableau)
                print(tableau)
                dual_simplex(tableau)
                flag=1
        print(tableau)'''

    def gomory_implementation(tableau,l):
        flag=0
        i=1
        while i<(tableau.shape[0]):
            f=tableau[i][0]%1
            #f=round(f,10)
            #print(f)
            #print('%.10f' % f)
            #if f>10^-7:
            #    print('True')
            if 1-10**(-7)>f and f>10**(-7):
                #print(l)
                (tableau,l)=gomory_tableau(tableau,l)
                #print(l)
                #print(tableau)
                (tableau,l)=dual_simplex(tableau,l)
                #print(tableau)
                #print()
                i=-1
            i=i+1
        #print(tableau)
        #print(l)
        return (tableau,l)

    '''def final_answer(tableau,l,n):
        #print(tableau)
        answer=[0]*n
        #print(l)
        for i in range (len(l)):
            if l[i]<=n:
                answer[l[i]-1]=round(tableau[i+1][0])
        return answer'''

    '''def gomory(filename):
        (n, m, c, A, b) = read_file(filename)
        #(table, basis, art) = convert_std(A, b, c, n, m)
        #print(table)
        (table,basis)=phase_2_simplex(A, b, c, n, m)
        #print(table)
        (table,basis)=gomory_implementation(table,basis)
        #print(table)
        ans=final_answer(table,basis,n)
        #print(ans)
        return ans'''

    A1=[[1,1],[-1,1],[6,2]]
    b1=[5,0,21]
    c1=[-2,-1]
    n1=2
    m1=3
    (n,m,c,A,b)=read_file(filename)
    #print(n,m,c,A,b)
    (table,basis)=phase_2_simplex(A, b, c, n, m)
    #print(f'table after 2 phase simplex:{table}')
    #print(table)
    #print(basis)
    (table,basis)=gomory_implementation(table,basis)
    #print(f'table after gomory:{table}')
    #print(f'basis after gomory:{basis}')
    #print(table)
    #print(basis)
    answer=[0]*n
        #print(l)
    """
    for i in range (len(basis)):
        if basis[i]<=n1:
            answer[basis[i]-1]=round(table[i+1][0])
    """
    for i in range(0, n):
        if (i+1) in basis:
            ind = basis.index(i+1)
            answer[i] = round(table[ind+1,0])
        else:
            answer[i] = 0
    return answer
    #print(final_answer(table,basis,n))
    # WE NOW TEST THE HELPER FUNCTIONS WE MADE FOR DRIVING OUT ARTIFICIAL VARIABLES AND FOR PHASE 2:
    """
    1. basic_art() - works
    2. check_zero() - works
    3. find_non_basic() - 
    4. find_next_row()
    5. remake_tableau() - ERROR
    6. re_reduce()
    7. phase_2_simplex()
    """
    """
    A=[[3,2],[-3,2]]
    b=[6,0]
    c=[0,-1]
    n=2
    m=2
    A1=[[1,1],[-1,1],[6,2]]
    b1=[5,0,21]
    c1=[-2,-1]
    n1=2
    m1=3
    A2=[[-1,1],[3,2],[2,3]]
    b2=[1,12,12]
    c2=[0,-1]
    n2=2
    m2=3
    """
    #(n, m, c, A, b) = read_file("input.txt")
    #print("printing answer")
    #print(phase_2_simplex(A1, b1, c1, n1, m1))
    #basis=[1,2,4]
    #tableau=np.array([[0.0,2,6,10,0,0],[2,-2,4,1,1,0],[-1,4,-2,-3,0,1]])
    #tableau1=np.array([[0.0,1,1,0,0],[-2,-1,-2,1,0],[-1,-1,0,0,1]])
    #dual_simplex(tableau1)
    #tableau2=np.array([[1.5,0,0,0.25,0.25],[1,1,0,1/6,-1/6],[1.5,0,1,0.25,0.25]])

    tableau3=np.array([[-31/4,0,0,-0.5,0,-0.25],[9/4,0,1,1.5,0,-0.25],[0.5,0,0,-2,1,0.5],[11/4,1,0,-0.5,0,0.25]])
    #tableau4=np.array([[-1.5,0,0,-0.25,-0.25],[1,1,0,1/6,-1/6],[1.5,0,1,0.25,0.25]])
    #gomory_tableau(tableau2)
    #(table,basis)=phase_2_simplex(A1, b1, c1, n1, m1)
    #print(basis)
    #print(table)
    #print(basis)

    #(table,basis)=gomory_implementation(table,basis)
    #print(table)
    #print(basis)
    #print(final_answer(table,basis,n))
    #print("printed answer")
    """
    (table, basis, art) = convert_std(A, b, c, n, m)
    print(art)
    print(table)

    print(is_finished(table))

    print(find_piv_row(table, 1))

    print(simplex_method(table, basis))
    print(basic_art(basis, art))
    print(check_zero(table, 1, art))
    """
#print(gomory("input.txt"))
