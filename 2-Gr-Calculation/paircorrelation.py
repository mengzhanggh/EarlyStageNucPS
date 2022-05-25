def pairCorrelationFunction_2D(x, y, S, rMax, dr):
    """Compute the two-dimensional pair correlation function, also known
    as the radial distribution function, for a set of circular particles
    contained in a square region of a plane.  This simple function finds
    reference particles such that a circle of radius rMax drawn around the
    particle will fit entirely within the square, eliminating the need to
    compensate for edge effects.  If no such particles exist, an error is
    returned. Try a smaller rMax...or write some code to handle edge effects! ;)

    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        S               length of each side of the square region of the plane
        rMax            outer diameter of largest annulus
        dr              increment for increasing radius of annulus

    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        annuli used to compute g(r)
        reference_indices   indices of reference particles
    """
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram
    # Number of particles in ring/area of ring/number of reference particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)

    # Find particles which are close enough to the box center that a circle of radius
    # rMax will not cross any edge of the box
    bools1 = x > rMax
    bools2 = x < (S - rMax)
    bools3 = y > rMax
    bools4 = y < (S - rMax)
    interior_indices, = where(bools1 * bools2 * bools3 * bools4)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a circle of radius rMax\
                will lie entirely within a square of side length S.  Decrease rMax\
                or increase the size of the square.")

    edges = arange(0., rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x) / S**2

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2)
        d[index] = 2 * rMax

        (result, bins) = histogram(d, bins=edges, normed=False)
        g[p, :] = result/numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (pi * (rOuter**2 - rInner**2))

    return (g_average, radii, interior_indices)
####

def pairCorrelationFunction_3D(x, y, z, S, rMax, dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell

    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        reference_indices   indices of reference particles
    """
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    bools1 = x > rMax
    bools2 = x < (S - rMax)
    bools3 = y > rMax
    bools4 = y < (S - rMax)
    bools5 = z > rMax
    bools6 = z < (S - rMax)

    interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
    #print(interior_indices)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x) / S**3

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2)
        d[index] = 2 * rMax

        (result, bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result / numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))

    return (g_average, radii, interior_indices)
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)
####


################################################################
################################################################
def pairCorrelationFunction_3D_meng(x, y, z, lowx,highx,lowy,highy,lowz,highz,deltaZ, rMax, dr, bd):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell

    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        reference_indices   indices of reference particles
        bd the bound of the box
    """
    import numpy as np
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    Sx=highx-lowx
    Sy=highy-lowy
    Sz=deltaZ
    #Sz=highz-lowz
    bools1 = x > bd[0]
    bools2 = x < bd[1]
    bools3 = y > bd[2]
    bools4 = y < bd[3]
    bools5 = z > bd[4]      
    bools6 = z < bd[5]

    interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
    #print(interior_indices)  #is a list of index of all particles
    num_interior_particles = len(interior_indices)
    print('total number of particle within the bound: ', len(x))
    print('number of particle within the bound: ', num_interior_particles)
    if num_interior_particles < 1:
        print("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")
        return ([0],[0],[0])

    edges = arange(0., rMax + 1.1 * dr, dr)
    #print(edges)
    ####    [0.  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.....................6 4.7 4.8 4.9 5.  5.1] shell
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments]) # empty  array for storing
    #print(g)
    #[[0. 0. 0. ... 0. 0. 0.] par1  num_interior_particles
    # [0. 0. 0. ... 0. 0. 0.] par2
    # [0. 0. 0. ... 0. 0. 0.].
    # ...
    # [0. 0. 0. ... 0. 0. 0.]
    # [0. 0. 0. ... 0. 0. 0.]
    # [0. 0. 0. ... 0. 0. 0.]]
    #   col  differet shells    num_increments
    radii = zeros(num_increments)
    numberDensity = len(x) / Sx*Sy*Sz   ### all particel / square box == mean density

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles): # loop over the interoir index
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2) # calculate the pair distance of current index to the rest of all points.            d is a list of distance

        d[index] = 2 * rMax ## the self distance is zero, to remove it, set it to 2 * rMax 

        (result, bins) = histogram(d, bins=edges, normed=False) ###  histogram gives the count in result, and the shells in bin

        #print('----')
        #print(histogram(d, bins=edges, normed=False))
        #  (array([ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        #        0,  0,  0,  6,  0,  0,  0,  0,  0,  0,  0, 12,  0,  0,  0,  0,  0,
        #        8,  0,  0,  0,  0,  0,  6,  0,  0,  0, 24,  0,  0,  0, 24,  0,  0]),            51 results
        #  array([0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. , 1.1, 1.2,
        #       1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. , 2.1, 2.2, 2.3, 2.4, 2.5,
        #       2.6, 2.7, 2.8, 2.9, 3. , 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8,
        #       3.9, 4. , 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5. , 5.1]))                 52 bins


        # base on the current particle z, calculate how much sphere volume outside the rectangle box, do adjustment

        adjustlist=meng_calc_outsideVol(bins, z[index], bd[4], bd[5])
        #print(adjustlist)
        #print(result)
        #print(len(adjustlist))
        #print(1-np.array(adjustlist))
        result_adj=result /(1-np.array(adjustlist))
        #print(result)
        g[p,:] = result_adj  # store the count for current shell current particle into g  row     
        #g[p,:] = result       # no adjsutment

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments) # empty  array for storing stat for each shell
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3) / numberDensity) 
######   result(count)/ (4/3 pi dr)    - current shell density
######--------------------
######   numberDensity                  - avg density 

    return (g_average, radii, interior_indices)
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)
####


def meng_calc_outsideVol(rlist,z, zlow,zhigh):
    if z < zlow or z > zhigh:
        print('z wrong --------------')
        print(z)
        print(zlow)
        print(zhigh)
        return
    else:
        outside_ratio=[]
        for i in range(len(rlist)-1): # rlist 52 element, range 0-50
            r = (rlist[i] + rlist[i+1]) / 2.
            height_top=(z+r)-zhigh
            outtop= calc_half_sphere(height_top ,r)
            #print(outtop)

            height_bot=zlow - (z-r)
            outbot= calc_half_sphere(height_bot ,r)
            #print(outbot)
            #print('----')
            outside_ratio.append(outtop + outbot)
            
        return outside_ratio


###https://keisan.casio.com/exec/system/1223382199 
def calc_half_sphere(h ,r):
    import math
    if h <=0:
        return 0
    else:
        dome_area=2 * math.pi * r * h
        sphera_total_area= 4 * math.pi * r**2
#        c=math.sqrt(h*(2*r - h))
#        v= math.pi / 6*h*(3*c**2 + h**2 )
        return dome_area/sphera_total_area










############################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
################################################################
def sphere(shape, radius, position):
    import numpy as np
    # assume shape and position are both a 3-tuple of int or float
    # the units are pixels / voxels (px for short)
    # radius is a int or float in px
    semisizes = (radius,) * 3
    #print(semisizes) ## (10, 10, 10)

    # genereate the grid for the support points
    # centered at the position indicated by position
    grid = [slice(-x0, dim - x0) for x0, dim in zip(position, shape)]
    position = np.ogrid[grid]

    #ogrid[0:5,0:5]
    #[array([[0],
    #        [1],
    #        [2],
    #        [3],
    #        [4]]), array([[0, 1, 2, 3, 4]])]


    # calculate the distance of all points from `position` center
    # scaled by the radius
    arr = np.zeros(shape, dtype=float)
    for x_i, semisize in zip(position, semisizes):
        # this can be generalized for exponent != 2
        # in which case `(x_i / semisize)`
        # would become `np.abs(x_i / semisize)`
        arr += (x_i / semisize) ** 2

    # the inner part of the sphere will have distance below 1
    return arr <= 1.0

def meng_calc_outsideVol2(xmin,xmax,ymin,ymax,zmin,zmax, rlist, ptx,pty,ptz): # rlist is a list of the shell radius
    import numpy as np
    # plot in 3D
    ## 100 z slice
    print('input parameter:')
    print(xmin,xmax,ymin,ymax,zmin,zmax, rlist, ptx,pty,ptz)

    outside_ratio=[]

    for i in range(len(rlist)-1): # rlist 52 element, range 0-50
        r1=rlist[i]
        r2=rlist[i+1]

        if r1==0:
            arr1=np.zeros((256,256,100), dtype=bool)
        else:
            arr1 = sphere((256,256, 100), r1, (ptx,pty,ptz)) ## ,    00011000
        arr2 = sphere((256,256, 100), r2, (ptx,pty,ptz)) ##          00111100


        arr=(np.invert(arr1==arr2))                               #  00100100   ## make the shell to be true in array
        shellvol= np.count_nonzero(arr)     #### the shell vol at specific radius
        if shellvol==0:
            outside_ratio.append(0)
            print('ratio:',0)
        else:
            ## put the rectangle cubic box in the volume
            arr[int(xmin):int(xmax),int(ymin):int(ymax),int(zmin):int(zmax)]=False
            shell_cut_vol= np.count_nonzero(arr)    ###  the shell vol at specific radius- substract the cubiod, which is the shell vol outside the box
            print('shell', i, '/',len(rlist)-1, '  ratioo:', shell_cut_vol /shellvol )  ### if no cubic cut, it will be equal to 1 
            outside_ratio.append(shell_cut_vol /shellvol)

    return outside_ratio

#            z,x,y = arr.nonzero()
#            import matplotlib.pyplot as plt
#            from mpl_toolkits.mplot3d import Axes3D
#            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d')
#            ax.scatter(x, y, z, c= 'red')

#    


def pairCorrelationFunction_3D_meng2(x, y, z, lowx,highx,lowy,highy,lowz,highz,deltaZ, rMax, dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell

    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        reference_indices   indices of reference particles
    """
    import numpy as np
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    ## high, low are the box boundary
    Sx=highx-lowx 
    Sy=highy-lowy
    Sz=deltaZ   ### deltaz were used to calculate the avg density
#    bools1 = x > lowx + rMax
#    bools2 = x < lowx + (Sx - rMax)
#    bools3 = y > lowy + rMax
#    bools4 = y < lowy + (Sy - rMax)
#    bools5 = z > lowz      
#    bools6 = z < highz

#    interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
    interior_indices, = where(x > 0) # all the x value must be >0, so it select all index
    print('hah')
    print(interior_indices)  #is a list of index of all particles    [   0    1    2 ... 1164 1165 1166]
   
    num_interior_particles = len(interior_indices)
    print(num_interior_particles)
    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax + 1.1 * dr, dr)
    #print(edges)
    ####    [0.  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.....................6 4.7 4.8 4.9 5.  5.1] shell
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments]) # empty  array for storing
    #print(g)
    #[[0. 0. 0. ... 0. 0. 0.] par1  num_interior_particles
    # [0. 0. 0. ... 0. 0. 0.] par2
    # [0. 0. 0. ... 0. 0. 0.].
    # ...
    # [0. 0. 0. ... 0. 0. 0.]
    # [0. 0. 0. ... 0. 0. 0.]
    # [0. 0. 0. ... 0. 0. 0.]]
    #   col  differet shells    num_increments
    radii = zeros(num_increments)
    numberDensity = len(x) / Sx*Sy*Sz   ### all particel / square box == mean density

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles): # loop over the interoir index
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2) # calculate the pair distance of current index to the rest of all points.   
        # d is a list of distance 

        d[index] = 2 * rMax ## the self distance is zero, to remove it, set it to 2 * rMax 

        (result, bins) = histogram(d, bins=edges, normed=False) ###  histogram gives the count in result, and the shells in bin
        #print(bins)
        #print('----')
        #print(histogram(d, bins=edges, normed=False))
        #  (array([ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        #        0,  0,  0,  6,  0,  0,  0,  0,  0,  0,  0, 12,  0,  0,  0,  0,  0,
        #        8,  0,  0,  0,  0,  0,  6,  0,  0,  0, 24,  0,  0,  0, 24,  0,  0]),            51 results, shell contained count
        #  array([0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. , 1.1, 1.2,
        #       1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. , 2.1, 2.2, 2.3, 2.4, 2.5,
        #       2.6, 2.7, 2.8, 2.9, 3. , 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8,
        #       3.9, 4. , 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5. , 5.1]))                 52 bins, shell radius


        # base on the current particle z, calculate how much sphere volume outside the rectangle box, do adjustment

        #adjustlist=meng_calc_outsideVol(bins, z[index], lowz, highz)
        print('------------particle#', p , '/',num_interior_particles)
        adjustlist=meng_calc_outsideVol2(lowx/11.68/2,
                                        highx/11.68/2,
                                        lowy/11.68/2,
                                        highy/11.68/2,
                                        lowz/11.68/2,
                                        highz/11.68/2,
                                        bins/11.68/2,
                                        x[index]/11.68/2,
                                        y[index]/11.68/2,
                                        z[index]/11.68/2)   # version two use the numeric np arry matrix volume calculation method

 


        #print(adjustlist)
        #print(result)
        #print(len(adjustlist))
        #print(1-np.array(adjustlist))
        result_adj=result /(1-np.array(adjustlist))  # 1-adjust means the ratio of shell vol inside the cubiod over the total shell vol
        #print(result)
        g[p,:] = result_adj  # store the count for current shell current particle into g  row     
        #g[p,:] = result       # no adjsutment

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments) # empty  array for storing stat for each shell
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3) / numberDensity) 
######   result(count)/ (4/3 pi dr)    - current shell density
######--------------------
######   numberDensity                  - avg density 

    return (g_average, radii, interior_indices)
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)
####



        

