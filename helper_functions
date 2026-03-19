###This file defines some helper functions
import matplotlib.pyplot as plt
import numpy as np
import sisl
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D


#Note -- the lines below are usefull to have in the jupyter notebook to force autoreload of the helper funcitons
#%load_ext autoreload
#%autoreload
#from helper_functions import *

#----------------------------------
def vdisplay(var):
    """Converts the var to a pretty string and inserts
    it on a new cell just below the present one.
    
    Then you have to change that 'next cell' type to Markdown and execute it.
    """
    # To print the var nicely.
    from pprint import pformat as pf

    string_to_insert=var
    
    # Create a code cell and insert a string in it
    get_ipython().set_next_input(string_to_insert)
    
    return
#----------------------------------


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)


def PrintMomentumRep(HamIn, use_tags=False,verbose=False,display_TeX=True,num_pf=False):
    """
    This function reads the hopping terms and prints momentum representation taht can be TeX printed
    """
    ###Extract all the hamiltonian elements
    HamElems = list(HamIn.iter_nnz())
    #print("HamElems:",HamElems)

    NO = HamIn.no ###Orbitals in unit cell
    ###Loop over elements in Hamiltonian and place them in the appropraite matrix place

    KElems={}
    ###Assign empty strings to all the indexes
    for x in range(NO):
        for y in range(NO):
            KElems[(x,y)]=''
    
    for x in HamElems:
        #print("-----")
        uc1 = np.mod(x[0],NO)
        uc2 = np.mod(x[1],NO)
        a1 = HamIn.o2a(uc1)
        a2 = HamIn.o2a(uc2)

        sc1 = HamIn.o2isc(x[0])
        sc2 = HamIn.o2isc(x[1])
        #print("x:",x)
        #print("uc1:",uc1)
        #print("uc2:",uc2)
        #print("sc1:",sc1)
        #print("sc2:",sc2)
        Tag1 = HamIn.geometry.atoms.orbital(x[0]).tag
        Tag2 = HamIn.geometry.atoms.orbital(x[1]).tag
        #print("Tag1:",Tag1)
        #print("Tag2:",Tag2)
        if not use_tags or Tag1 == '':
            Tag1 = uc1
        if not use_tags or Tag2 == '':
            Tag2 = uc2

        Rdiff = sc2-sc1
        #print("Rdiff:",Rdiff)
        HamElemAB = HamIn[x[0],x[1]]
        #print("HamElemAB:",HamElemAB)
        #print(PreFactor)
        if np.sum(np.abs(Rdiff)) != 0:
            #print("Add momentum part")
            Rvalue=''
            MomPart="e^{"
            for dim in range(3):
                Rvec="\mathbf{{R}}_{dim}".format(dim=dim+1)
                kvec="\mathbf{{k}}_{dim}".format(dim=dim+1)
                Rkvec=Rvec+kvec
                if Rdiff[dim] == 1:
                    MomPart += "-\imath {Rkvec}".format(Rkvec=Rkvec)
                    Rvalue += "+"+Rvec
                elif Rdiff[dim] == -1:
                    MomPart += "+\imath {Rkvec}".format(Rkvec=Rkvec,sc=-Rdiff[dim])
                    Rvalue += "-"+Rvec
                elif Rdiff[dim] > 0:
                    MomPart += "-\imath {sc} {Rkvec}".format(Rkvec=Rkvec,sc=-Rdiff[dim])
                    Rvalue += "+ {sc} {Rvec}".format(Rvec=Rvec,sc=Rdiff[dim])
                    
                elif Rdiff[dim] < 0:
                    MomPart += "+\imath {sc} {Rkvec}".format(Rkvec=Rkvec,sc=Rdiff[dim])
                    Rvalue += "- {sc} {Rvec}".format(Rvec=Rvec,sc=-Rdiff[dim])
            MomPart +="}"
            #print(MomPart)
        else:
            MomPart=""
            Rvalue='0'
        if num_pf:
            PreFactor=str(HamElemAB)
        else:
            PreFactor="a_{{{uc1},{uc2}}}^{{({Rvalue})}}".format(uc1=Tag1,uc2=Tag2,Rvalue=Rvalue)

            
        Factor=PreFactor+MomPart
        #print(Factor)
        if KElems[(uc1,uc2)] == '':
            KElems[(uc1,uc2)] = Factor
        else:
            KElems[(uc1,uc2)]=KElems[(uc1,uc2)]+"+"+Factor

    #print(KElems)
    ###Now we build the actual matrix
    TeXString="\\begin{pmatrix}"
    for x in range(NO):
        for y in range(NO):
            if KElems[(x,y)] != "":
                TeXString+=KElems[(x,y)]
            else:
                TeXString+="0"
            if y+1 != NO:
                TeXString+=" & "
            else:
                if x+1 != NO:
                    TeXString+=" \\\\ " ##Add a newline
    TeXString+="\\end{pmatrix}"
    #print(TeXString)
    #vdisplay(TeXString)
    #return TeXString
    ####Display the text as markdown!
    from IPython.display import display, Markdown
    if display_TeX:
        display(Markdown(TeXString))            
    return TeXString


def draw_unit_cell(Lattice,fig,ucc=[0,0,0]):
    Vec1 = Lattice.vertices()[1,0,0]
    Vec2 = Lattice.vertices()[0,1,0]
    Vec3 = Lattice.vertices()[0,0,1] 

    V12 = Vec1+Vec2
    Shift = -V12/2+ucc

    ax = fig.gca()
    ax.plot([Shift[0],Shift[0]+Vec1[0]],[Shift[1],Shift[1]+Vec1[1]],"--r")
    ax.plot([Shift[0],Shift[0]+Vec2[0]],[Shift[1],Shift[1]+Vec2[1]],"--r")
    ax.plot([Shift[0]+Vec1[0],Shift[0]+Vec1[0]+Vec2[0]],
             [Shift[1]+Vec1[1],Shift[1]+Vec1[1]+Vec2[1]],"--r")
    ax.plot([Shift[0]+Vec2[0],Shift[0]+Vec1[0]+Vec2[0]],
             [Shift[1]+Vec2[1],Shift[1]+Vec1[1]+Vec2[1]],"--r")


def draw_atoms(SISGeometry,ax,AtomRadius=None,zorder=2,dim='2D',color='red'):
    NA_S = SISGeometry.na_s ####Number of atoms in super cell
    if AtomRadius is None:
        MinDist = GetMinDistance(SISGeometry)
        AtomRadius = MinDist*0.2
    for a_no in range(NA_S):
        #print("Atom number",a_no) if verbose else False
        Position =SISGeometry.axyz(a_no) ### The super cell offset
        #print("xyz of elem:",Position[0],Position[1]) if verbose else False

        color_atom=get_cyclic_option(color,a_no)
        radius_atom=get_cyclic_option(AtomRadius,a_no)

        
        if dim == '2D':
            circle = plt.Circle((Position[0],Position[1]), radius_atom, color=color_atom,zorder=zorder)
            ax.add_patch(circle)
            ###This line ensures that the x/y limits are aware of the positons
            ###The circle patch is not nececarrily picked up by the plot
            ax.plot(Position[0],Position[1],"-")

        elif dim == '3D':
            sphere = ax.scatter(Position[0],Position[1],Position[2],s=100*AtomRadius**2,color=color_atom)
 

        

def PlotHamTerms(HamIn, AtomRadius=None,uc=False, Tags=None,Atoms=None,DrawType=None,
                 cmap="hsv",ucc=[0,0,0],hold_show=False,verbose=False,dim='2D',
                 hw=2,hl=5,tw=1.0,color='red',figsize=[8,5]):
    """
    Plots the terms of a Hamiltonian in a 2D representation using the x and y components.

    Args:
        HamIn (object):
                       Input Hamiltonian object containing geometric and Hamiltonian data.
        uc (bool, optional):
                       If True, draws the unit cell in the plot. Defaults to False.
        Tags (list, optional):
                       List of tags to filter which atoms to show. If None, shows all atoms. Defaults to None.
        Atoms (list, optional):
                       List of atom indices to include in the plot. If None, includes all atoms. Defaults to None.
        AtomRadius (float, optional):
                       Radius for the visual representation of atoms. Defaults to None.
        DrawType (str, optional):
                       Specifies how to draw the Hamiltonian terms. Options are 'Phase' or other (default will represent magnitudes).
        cmap (str, optional):
                       Colormap to use for visual representation of phases. Defaults to "hsv".
        ucc (list, optional):
                       Coordinates for the unit cell center. Defaults to [0, 0, 0].
        hold_show (bool, optional):
                       If False, displays the plot immediately. If True, holds the display. Defaults to False.
        verbose (bool, optional):
                       If True, prints additional debugging information. Defaults to False.

    Returns:
        matplotlib.figure.Figure: The matplotlib figure object with the plotted Hamiltonian terms.

    Notes:
        - The function iterates over all non-zero elements of the Hamiltonian and draws arrows representing the interactions based on the Hamiltonian terms.
        - The drawn arrows can represent the phase or magnitude of the Hamiltonian terms depending on the `DrawType` specified.

    Example:
        fig = PlotHamTerms(my_hamiltonian, verbose=True, Tags=["A", "B"], DrawType="Phase")
    """
    
    ##print('hello')
    fig = plt.figure(figsize=figsize)
    if dim=='2D':
        ax = fig.gca()
    elif dim=='3D':
        ax = fig.add_subplot(projection='3d')
    else:
        raise ValueError("'dim' needs to be '2D' or '3D'")

    ###Loop over atoms in extened unit cell and draw them

    ####Get the Ham elements
    HamElems = list(HamIn.iter_nnz())
    print(HamElems)  if verbose else False

    ###Draw the atoms
    draw_atoms(HamIn.geometry,ax,AtomRadius=AtomRadius,zorder=2,dim=dim,color=color)

    AtomsPerCell = HamIn.na
    OrbtalsPerCell = HamIn.no
    ###Loop over elements in Hamiltonian
    for x in HamElems:
        if x[0] == x[1]:
            ###If the two orbitals are the same. Skip it.
            continue
        if Tags != None: ###Only show atoms on the tags list
            Tag1 = HamIn.geometry.atoms.orbital(x[0]).tag
            Tag2 = HamIn.geometry.atoms.orbital(x[1]).tag
            if Tag1 not in Tags:
                continue
            if Tag2 not in Tags:
                continue
        print("x of elem:",x) if verbose else False
        a1 = HamIn.o2a(x[0])
        a2 = HamIn.o2a(x[1])
        if Atoms != None: ###Only show atoms on the atom list
            a1uc = np.mod(a1,AtomsPerCell)
            a2uc = np.mod(a2,AtomsPerCell)
            if a1uc not in Atoms:
                continue
            if a2uc not in Atoms:
                continue
        print("  a1:",a1)  if verbose else False
        print("  a2:",a2)  if verbose else False
        r1 = HamIn.axyz(a1)
        r2 = HamIn.axyz(a2)
        #plt.plot([r1[0],r2[0]],[r1[1],r2[1]],"-r")
        Differnce=r2-r1
        #print("Differnce:", Differnce) if verbose else False

        HamElem = HamIn[x[0],x[1]]
        
        print("x[0]",x[0])   if verbose else False
        print("x[1]",x[1])   if verbose else False
        print("HamElem:",HamElem)   if verbose else False

        ###FIXME this is a workaround for SOC where the Hamiltonian has 8 components!
        try: ###Lists and scalar need to be handles differnes
            SumTerms = sum(abs(HamElem))
        except TypeError as te:
            SumTerms = HamElem
       
        if SumTerms == 0.0:
            continue
        
        

        ###Get the colorbal for bwr
        cmapSample = plt.get_cmap(cmap)
        
        ###Loop over the super cells and find the start atom
        for sc_indx in range(HamIn.n_s):
            print("  sc_indx:",sc_indx)  if verbose else False
            a1sc = a1+AtomsPerCell*sc_indx
            print("  a1sc:",a1sc)  if verbose else False
            if (a1sc < HamIn.na_s):
                r1sc = HamIn.axyz(a1sc)
                r2sc = r1sc + Differnce
                #print(r1,r2)  if verbose else False

                if DrawType == "Phase":
                    Phase = np.mod(np.angle(HamElem),2*np.pi)/(2*np.pi) ####Number in range 0 - 1
                    PhaseColor = cmapSample(Phase)
                    plt.arrow(r1sc[0],r1sc[1],Differnce[0],Differnce[1],width=0.1,color=PhaseColor,
                              length_includes_head=True,zorder=3)
                else:
                    r1sc2D=1*r1sc;r2sc2D=1*r2sc
                    r1sc2D[2]=0;r2sc2D[2]=0
                    ArrLen=np.linalg.norm(r2sc2D-r1sc2D)
                    #plt.arrow(r1sc[0],r1sc[1],Differnce[0],Differnce[1],
                    #          width=0.03*ArrLen,
                    #          length_includes_head=True,zorder=3)
                    #plt.plot([r1[0],r2[0]],[r1[1],r2[1]],"-r")


                    style = ("Simple, tail_width="+str(tw)+", "+
                             "head_width="+str(hw)+", "+
                             "head_length="+str(hl))
                    kw = dict(arrowstyle=style, color="k")
                    if dim=='2D':
                        a3 = patches.FancyArrowPatch((r1sc[0],r1sc[1]),
                                                     (r2sc[0],r2sc[1]),
                                                     connectionstyle="arc3,rad=.2",
                                                     **kw,zorder=3)
                        plt.gca().add_patch(a3)
                    elif dim=='3D':
                        a3 = Arrow3D((r1sc[0],r2sc[0]),
                                     (r1sc[1],r2sc[1]),
                                     (r1sc[2],r2sc[2]),
                                     connectionstyle="arc3,rad=.2",
                                     **kw)
                        plt.gca().add_artist(a3)
                    
                        

    if uc:
        draw_unit_cell(HamIn.lattice,fig,ucc=ucc)

    plt.xlabel("x-position [Å]")
    plt.ylabel("y-position [Å]")
        
    if dim=='2D':
        plt.axis("equal")
    elif dim=='3D':
        ax.set_zlabel("z-position [Å]")
        plt.axis("auto")
    if not hold_show:
        plt.show(block=False)
    return fig

    ######### End of function




def CheckHermiticity(HamIn, verbose=False,succes_message="Hamiltonian is hermitian!"):
    """
    This function checks that the hamiltonian is Hermitian!
    """

    ###Extract all the hamiltonian elements
    HamElems = list(HamIn.iter_nnz())
    #print("HamElems:",HamElems)

    NO = HamIn.no ###Orbitals in unit cell
    ###Loop over elements in Hamiltonian and check that there is acorresponding term back as well
    for x in HamElems:
        #print("-----")
        fuc1 = np.mod(x[0],NO)
        fuc2 = np.mod(x[1],NO)
        a1 = HamIn.o2a(fuc1)
        a2 = HamIn.o2a(fuc2)
        fsc1 = HamIn.o2isc(x[0])
        fsc2 = HamIn.o2isc(x[1])
        #print("x:",x)
        #print("tuc1:",fuc1)
        #print("tuc2:",fuc2)
        #print("tsc1:",fsc1)
        #print("tsc2:",fsc2)
        fdfs = fsc2-fsc1
        #print("fdfs:",fdfs)
        HamElemAB = HamIn[x[0],x[1]]
        #print("HamElemAB:",HamElemAB)
        FoundMatch=False
        for y in HamElems:
            #print("  -.-.-.-.-")
            tuc1 = np.mod(y[0],NO)
            tuc2 = np.mod(y[1],NO)
            tsc1= HamIn.o2isc(y[0])
            tsc2= HamIn.o2isc(y[1])
            #print("  y:",y)
            #print("  tuc1:",tuc1)
            #print("  tuc2:",tuc2)
            #print("  tsc1:",tsc1)
            #print("  tsc2:",tsc2)
            tdfs = tsc1-tsc2 ###Reverse diff
            #print("  tdfs:",tdfs)
            HamElemBA = HamIn[y[0],y[1]]
            #print("  HamElemBA:",HamElemBA)

            if fuc1 != tuc2:
                continue
            if fuc2 != tuc1:
                continue
            #print("Same orbital hop -!")
            if np.sum(np.abs(fdfs - tdfs)):

                ###If this is not zero they are not the same
                continue
            #print("Same super cell differnces, now we can compare!")
            if np.isclose(HamElemAB,np.conj(HamElemBA)):
                ###We found a match, we are done!
                #print("We found a match. Let's break")
                FoundMatch=True
                break
            else:
                print(HamElemAB,HamElemBA)
                StringText = 'Hamiltonian is not Hermitian! Orbital {OA} (on atom {aA} ) and {OB} (on atom {aB})  have hopping  terms {OA} -> {OB} {ShiftAB} = {HAB} and {OB} -> {OA} {ShiftBA} = {HBA} respectively!'
                StringError=StringText.format(OA=fuc1,OB=fuc2,ShiftAB=fdfs,ShiftBA=-fdfs,
                                              HAB=HamElemAB,HBA=HamElemBA,aA=a1,aB=a2)
                raise ValueError(StringError)
        ### If we to not find a match that is also and error
        if not FoundMatch:
            #print("Ops we did not!")
            StringText = 'Hamiltonian is not Hermitian! There is a term missing for Orbital {OA} (on atom {aA} ) to {OB} (on atom {aB}) and sc shift {ShiftAB}. We have hopping terms {OB} -> {OA} {ShiftAB} = {HAB} but not the one back: {OB} -> {OA} {ShiftBA}!'
            StringError=StringText.format(OA=fuc1,OB=fuc2,ShiftAB=fdfs,
                                          ShiftBA=-fdfs,
                                          HAB=HamElemAB,aA=a1,aB=a2)
            print(StringError)
            raise ValueError(StringError)
    if succes_message != None:
        print(succes_message)

    
def GetMinDistance(GeoIn):
    NA_S = GeoIn.na_s ####Number of atoms in super cell
    AtomPos=np.zeros((NA_S,3))
    for a_no in range(NA_S):
        Position=GeoIn.axyz(a_no) ### The super cell offset
        AtomPos[a_no,:] = Position[:3]
    #display(AtomPos)
    #print("MinDist:",MinDist)
    Vec1 = GeoIn.lattice.vertices()[1,0,0]
    Vec2 = GeoIn.lattice.vertices()[0,1,0]
    Vec3 = GeoIn.lattice.vertices()[0,0,1]
    A1=np.linalg.norm(Vec1)
    A2=np.linalg.norm(Vec2)
    A3=np.linalg.norm(Vec3)

    MinDist=np.min([A1,A2,A3])
    #print("Min Distance:",MinDist)
    
    for n1 in range(NA_S):
        for n2 in range((n1+1),NA_S):
            #print("------------")
            #print("n1,n2:",n1,n2)
            A1=AtomPos[n1,:]
            A2=AtomPos[n2,:]
            #print("A1,A2:",A1,A2)
            Dist=np.linalg.norm(A2-A1)
            #print("Dist:",Dist)
            if Dist > 0:
                if MinDist == None or Dist < MinDist:
                    MinDist = Dist
    return MinDist

    
def PlotGeometry(GeoIn, verbose=False,ucc=[0,0,0],uc=False,hold_show=False,AtomRadius=None,dim='2D',color='red',figsize=[8,5]):
    fig = plt.figure(figsize=figsize)
    ###Loop over atoms in extened unit cell and draw them

    NA_S = GeoIn.na_s ####Number of atoms in super cell
    print("Atoms in super cell",NA_S) if verbose else False
    XYZ = GeoIn.xyz ###Position of atoms within the cell
    print("Positions in cell:\n",XYZ) if verbose else False

    if dim=='2D':
        ax = fig.gca()
    elif dim=='3D':
        ax = fig.add_subplot(projection='3d')
    else:
        raise ValueError("'dim' needs to be '2D' or '3D'")
    ###Draw the atoms
    draw_atoms(GeoIn,ax,AtomRadius=AtomRadius,zorder=2,dim=dim,color=color)

    if uc:
        draw_unit_cell(GeoIn.lattice,fig,ucc=ucc)

        
    plt.xlabel("x-position [Å]")
    plt.ylabel("y-position [Å]")
    if dim=='2D':
        plt.axis("equal")
    elif dim=='3D':
        ax.set_zlabel("z-position [Å]")
        plt.axis("auto")

    if not hold_show:
        plt.show(block=False)
    return fig, ax

    ######### End of function

def EnergySurface(MPG,EigsMPG=None,verbose=False,Ham=None,figsize=[8,5]):
    """
    Make a surface plot from a given Monkhorst grid of the brilluin zone
    MPG = sisl.physics.MonkhorstPack(H, [Nx, Ny, 1])
    """

    from matplotlib import cm
    from matplotlib.ticker import LinearLocator

    if EigsMPG is None: ###Compute the eigenvalues for all the k-points
        print("Eigenvalues not known ... computing them ...")
        EigsMPG  = MPG.apply.array.eigh()
        print("Should you already have them precomuted, then you can supply them, with the EigsMPG=MPG.apply.array.eigh()")

        
    kxlist=MPG.k[:,0]
    kylist=MPG.k[:,1]
    
    kxUnique=np.sort(np.unique(kxlist))
    kyUnique=np.sort(np.unique(kylist))

    print(kxUnique) if verbose else False
    print(kyUnique) if verbose else False
  
    NoSurfaces=len(EigsMPG[0])
    #print("NoSurfaces:",NoSurfaces) #if verbose else False
    
    Egrid = np.zeros([len(kxUnique),len(kyUnique),NoSurfaces])*np.nan

    for n in range(len(kxlist)):
        ####slow but workds
        xpos = np.argmax(kxUnique == kxlist[n])
        ypos = np.argmax(kyUnique == kylist[n])
        ##print(xpos,ypos)
        Egrid[xpos,ypos,:]=EigsMPG[n]



    vmin=np.min(Egrid)
    vmax=np.max(Egrid)
    nk1, nk2 = np.meshgrid(kxUnique, kyUnique)


    if Ham is not None:
        #print("Reading reciprocal lattice from Ham")
        ###Extract the reciprocal lattice vectors
        K1 = Ham.lattice.icell[0,:]
        K2 = Ham.lattice.icell[1,:]
        K3 = Ham.lattice.icell[2,:]

        #### NewK =  K1*nk1 + K2*nk2 ....

        ### X -coords = K1[0]*nk1 + K2[0]*nk2 + K3[0]*nk3
        ### Y -coords = K1[1]*nk1 + K2[1]*nk2 + K3[1]*nk3
        ### Z -coords = K1[2]*nk1 + K2[2]*nk2 + K3[2]*nk3

        X = K1[0]*nk1 + K2[0]*nk2 ###No nk3
        Y = K1[1]*nk1 + K2[1]*nk2 ###No nk3
    else:
        X = nk1
        Y = nk2

    
    # Plot the surface.
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=figsize)
    for surfno in range(NoSurfaces):
        surf = ax.plot_surface(X, Y, Egrid[:,:,surfno], cmap=cm.coolwarm,
                               linewidth=0, antialiased=False,vmin=vmin,vmax=vmax)
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    if Ham is not None:
        plt.xlabel("kx [1/Å]")
        plt.ylabel("ky [1/Å]")
    else:
        plt.xlabel("nkx - coord allong K_1")
        plt.ylabel("nky - coord allong K_2")
    ax.set_zlabel("Energy [eV]")

    plt.show(block=False)


    dx=np.diff(kxUnique)[0]
    dy=np.diff(kyUnique)[0]
    ###Use meshgrid
    ##XM,YM = np.meshgrid(kxUnique,kyUnique)
    #print(XM)

    Extent=[np.min(kxUnique)-dx/2,np.max(kxUnique)+dx/2,
            np.min(kyUnique)-dy/2,np.max(kyUnique)+dy/2]

        
    if False:
        plt.figure()
        plt.imshow(np.flipud(np.transpose(Egrid[:,:,0])),extent=Extent)
        plt.xlabel("kx - in units of eigenvector")
        plt.ylabel("ky - in units of eigenvector")
        plt.axis("equal")
        plt.colorbar()
        plt.show(block=False)
        
        #input("Press enter to stop")
    return fig, ax
    


####FIXME --- Med the DOS and option, not a standard

def plot_bands(kpath, ymin=None, ymax=None,Esteps=101,bandcolor=None,show_band=None,linestyle='solid',linewidth=1,hold_show=False,DOSKpoints=[100,100,1],show_dos=False,H=None,figure=None):
    """
    Generate a band structure and density of states (DOS) plot from k-path data.

    Parameters:
    ----------
    kpath : object
        An object containing information about the k-path, including methods
        to retrieve eigenvalues and ticks for plotting.
    
    ymin : float, optional
        Minimum value for the energy axis; if not provided, defaults to the 
        lower limit of the eigenvalues.
    bandcolor : string or list of strings
        Choose your colors of the bands
    ymax : float, optional
        Maximum value for the energy axis; if not provided, defaults to the 
        upper limit of the eigenvalues.

    Returns:
    -------
    None
        Displays a figure with two subplots: one for the band structure 
        and one for the DOS plot.
    """
    
    eigs = kpath.apply.array.eigh()
    ###Select only some bands
    if show_band != None:
        eigs = eigs[:,show_band]

    if show_dos:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize =(11, 5))
        ax = ax1, ax2
    else:
        if figure == None:
            fig, (ax1) = plt.subplots(1, 1, figsize =(8, 5))
        else:
            fig = figure
            ax1 = fig.gca()
        ax = [ax1]
        
        

    # Retrieve the tick-marks and the linear k points
    xtick, xtick_label = kpath.lineartick()
    lk = kpath.lineark()
    for band_id in range(eigs.shape[1]):
        color=get_cyclic_option(bandcolor,band_id)
        lineStyle=get_cyclic_option(linestyle,band_id)
        lineWidth=get_cyclic_option(linewidth,band_id)
            
        ax1.plot(lk, eigs[:,band_id],color=color,linestyle=lineStyle,linewidth=lineWidth) # Give each band a given

    ax1.set_xticks(xtick)
    ax1.set_xticklabels(xtick_label)
    ax1.set_xlim(0, lk[-1])
    YMIN,YMAX=ax1.get_ylim()
    if ymin != None:
        YMIN=ymin
    if ymax != None:
        YMAX=ymax
    ax1.set_ylim(YMIN,YMAX)

    En = np.linspace(YMIN, YMAX, Esteps) # Energy array to calculate the DOS from
    NumBands=len(eigs[0])

    if show_dos:
        if H ==  None:
            raise ValueError("Please supply the hamiltonians with plot_bands(kpath,H=<Ham>,...)")
            
        
        DOSTotal=np.zeros(len(En))
        ##FIXME sometimes the gaussian is not wide enough... we can mofiy this using the distribution
        ###https://sisl.readthedocs.io/en/v0.14.3/api/generated/sisl.physics.get_distribution.html

        ###Compute the DOS from all points in BZ
        bz = sisl.physics.MonkhorstPack(H, DOSKpoints) # We sample the BZ
        MHEigs = bz.apply.array.eigh()
        for bandno in range(NumBands):
            DOSforband=sisl.physics.electron.DOS(En,MHEigs[:,bandno])
            #ax2.plot(DOSforband, En) # Plot the DOS FIX MEEEEE
            DOSTotal += DOSforband

        ax2.plot(DOSTotal, En) # Plot the DOS FIX 
        ax2.set_xlabel('DOS [Arb units]')

        XMIN,XMAX=ax2.get_xlim()
        ax2.set_xlim(0,XMAX)
        ax2.set_xticks([0])

        
    # Also plot x-major lines at the ticks
    for tick in xtick:
        ax1.plot([tick,tick], [YMIN, YMAX], 'k')


    
    ax1.set_xlabel('Momentum')

    for i in ax:
        i.set_ylim(YMIN, YMAX)
        i.set_ylabel('Energy [eV]')

    if not hold_show:
        plt.show(block=False)    
    return fig,ax ###return the future future handling


def plot_spin_bands(kpath, bands, spin_moments,ymin=None,ymax=None,
                    figsize=(8, 4.5)):

    from matplotlib.collections import LineCollection
    
    lk = kpath.lineark()
    xtick, xtick_label = kpath.lineartick()    
    nk, nbands = bands.shape

    YMIN=np.min(bands)
    YMAX=np.max(bands)

    
    # Create a figure with three subplot one for each component of the spin moment
    fig, axes = plt.subplots(1, 3, figsize=figsize, dpi=300, sharex=True, sharey=True)

    # Set the range of z-values, which will determine the color.
    norm = plt.Normalize(-1, 1)

    # Iterate over the spin components
    for icomp, component in enumerate(['$S_x$', '$S_y$', '$S_z$']):
        # Iterate over all bands 
        for ibnd in range(nbands):
            # It is not possible to change the color of a line directly, so we create small 
            # line segements from one point on the x-axis to the next. These segments can 
            # then be colored individually.
            points = np.array([lk, bands[:, ibnd]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            # Create a collection of the segments and specify a map that assigns colors 
            # to the segments according to the z-value
            lc = LineCollection(segments, cmap='coolwarm', norm=norm)
            # Set the z-values 
            lc.set_array(spin_moments[:, ibnd, icomp])             
            lc.set_linewidth(3)

            # Add the LineCollection to the subplot
            line = axes[icomp].add_collection(lc)
            
        axes[icomp].set_title(component)
        axes[icomp].set_xlabel('Momentum')


    # All subplots share the same axis settings, so we can just adjust them once

    if ymin != None:
        YMIN=ymin
    if ymax != None:
        YMAX=ymax
    axes[0].set_ylim(YMIN,YMAX)

    axes[0].set_xlim(min(lk), max(lk))
    axes[0].set_ylabel('Eigenspectrum [eV]')
    axes[0].xaxis.set_ticks(xtick)
    axes[0].set_xticklabels(xtick_label)

    for axis in axes:
        for tick in xtick:
            axis.plot([tick, tick], [ymin, ymax], 'k', linewidth=0.5)
    
    # Add a colorbar to the plot
    fig.colorbar(line, ax=axes.ravel().tolist())
    plt.show(block=False)
    return fig, axes



def get_cyclic_option(DataList,indx):
    if type(DataList) == list or type(DataList) == np.ndarray:
        indx_id = np.mod(indx,len(DataList))
        return DataList[indx_id]
    else:
        return DataList
