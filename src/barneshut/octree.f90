! ************************************************************
!! Octree (Barnes-Hut)
!
!> Objectives
!  This module generates a octree based in the positions of the
!  particles and evaluates the forces using the Barnes-Hut 
!  approximation with a quadrupole expansion.
!
!> Modified
!  2026.06.15
!
!> Author
!  oap
!
MODULE octree_mod
    USE types_mod

    IMPLICIT NONE
    PRIVATE
    PUBLIC OctreeType

    TYPE :: OctreeType
        ! max depth of the tree
        INTEGER :: max_depth = 32
        ! amplificator of the size of the quadtree-root
        REAL(pf) :: side_amplificator = 1.2_pf

        ! to save_txt
        INTEGER :: save_txt

        ! masses, positions and number of bodies (N)
        REAL(pf), ALLOCATABLE :: m(:), x(:), y(:), z(:)
        INTEGER :: N

        ! nodes
        INTEGER :: number_of_nodes = 0
        INTEGER :: max_number_of_nodes
        REAL(pf), ALLOCATABLE :: ns_cx(:), ns_cy(:), ns_cz(:) ! centers
        REAL(pf), ALLOCATABLE :: ns_halfside(:), ns_L2(:)
        REAL(pf), ALLOCATABLE :: ns_mass(:), ns_qcm_x(:), ns_qcm_y(:), ns_qcm_z(:)
        REAL(pf), ALLOCATABLE :: ns_quad(:,:)
        INTEGER, ALLOCATABLE :: ns_particle(:)
        INTEGER, ALLOCATABLE :: ns_type(:)
        INTEGER, ALLOCATABLE :: ns_child(:,:)
        INTEGER, ALLOCATABLE :: ns_depth(:)
    CONTAINS
        PROCEDURE :: init
        PROCEDURE :: allocate_nodes, add_node, allocate_subnode, add_to_subnode, add
        PROCEDURE :: forces => evaluate_forces_over_p
        PROCEDURE :: evaluate_quad, evaluate_quad_aprox
    END TYPE

CONTAINS

SUBROUTINE init (self, m, x, y, z, save_txt)
! this subroutine inits the tree by allocating the global vectors and adding each particle
! in a node. if its the case it saves the root information too.
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: m(:), x(:), y(:), z(:)
    INTEGER, OPTIONAL :: save_txt
    REAL(pf) :: infos_root(4)
    INTEGER :: p, idx_root

    ! saving particles information
    self % N = SIZE(m)
    ALLOCATE(self % m(self % N))
    ALLOCATE(self % x(self % N))
    ALLOCATE(self % y(self % N))
    ALLOCATE(self % z(self % N))
    self % m = m
    self % x = x
    self % y = y
    self % z = z

    ! init the root
    self % max_number_of_nodes = 8 * self % N
    self % number_of_nodes = 0
    CALL self % allocate_nodes()
    
    infos_root = node_size_center(x, y, z)
    CALL self % add_node(infos_root(1), infos_root(2), infos_root(3), 1.2*infos_root(4), 0, idx_root)

    ! if wants to save_txt
    self % save_txt = -1
    IF (PRESENT(save_txt)) THEN
        self % save_txt = save_txt
        WRITE (self % save_txt, *) self%ns_cx(1), self%ns_cy(1), self%ns_cz(1), self%ns_halfside(1)
    ENDIF

    ! add the children to root
    DO p = 1, self % N
        CALL self % add(idx_root, p)
    END DO
END SUBROUTINE

SUBROUTINE allocate_nodes (self)
! this subroutine allocates the global vectors that stores informations about the nodes.
! if the tree already have nodes, so its a reallocation. in that case, the old information
! are stored in temp vectors and the vectors are reallocated with twice times the old size
! without lose the old information.
! @calledby init, add_node
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), ALLOCATABLE :: temp_real(:)
    INTEGER, ALLOCATABLE :: temp_int(:), temp_int_2(:,:)
    INTEGER :: old_size, new_size

    ! if the tree already exists, so is the case of reallocation
    IF (self % number_of_nodes > 0) THEN
        old_size = self % max_number_of_nodes
        self % max_number_of_nodes = 2 * old_size
        new_size = self % max_number_of_nodes
        
        ! allocate the temp vectors
        ALLOCATE(temp_real(old_size))
        ALLOCATE(temp_int(old_size))
        ALLOCATE(temp_int_2(old_size,8))

        ! now deallocate and reallocate
        temp_real = self % ns_cx
        DEALLOCATE(self % ns_cx)
        ALLOCATE(self % ns_cx(new_size))
        self % ns_cx(1:old_size) = temp_real

        temp_real = self % ns_cy
        DEALLOCATE(self % ns_cy)
        ALLOCATE(self % ns_cy(new_size))
        self % ns_cy(1:old_size) = temp_real

        temp_real = self % ns_cz
        DEALLOCATE(self % ns_cz)
        ALLOCATE(self % ns_cz(new_size))
        self % ns_cz(1:old_size) = temp_real

        temp_real = self % ns_halfside
        DEALLOCATE(self % ns_halfside)
        ALLOCATE(self % ns_halfside(new_size))
        self % ns_halfside(1:old_size) = temp_real

        temp_real = self % ns_L2
        DEALLOCATE(self % ns_L2)
        ALLOCATE(self % ns_L2(new_size))
        self % ns_L2(1:old_size) = temp_real

        temp_real = self % ns_mass
        DEALLOCATE(self % ns_mass)
        ALLOCATE(self % ns_mass(new_size))
        self % ns_mass(1:old_size) = temp_real

        temp_real = self % ns_qcm_x
        DEALLOCATE(self % ns_qcm_x)
        ALLOCATE(self % ns_qcm_x(new_size))
        self % ns_qcm_x(1:old_size) = temp_real

        temp_real = self % ns_qcm_y
        DEALLOCATE(self % ns_qcm_y)
        ALLOCATE(self % ns_qcm_y(new_size))
        self % ns_qcm_y(1:old_size) = temp_real

        temp_real = self % ns_qcm_z
        DEALLOCATE(self % ns_qcm_z)
        ALLOCATE(self % ns_qcm_z(new_size))
        self % ns_qcm_z(1:old_size) = temp_real

        temp_int = self % ns_particle
        DEALLOCATE(self % ns_particle)
        ALLOCATE(self % ns_particle(new_size))
        self % ns_particle(1:old_size) = temp_int

        temp_int = self % ns_type
        DEALLOCATE(self % ns_type)
        ALLOCATE(self % ns_type(new_size))
        self % ns_type(1:old_size) = temp_int
        
        temp_int_2 = self % ns_child
        DEALLOCATE(self % ns_child)
        ALLOCATE(self % ns_child(new_size, 8))
        self % ns_child(1:old_size,:) = temp_int_2

        DEALLOCATE(temp_real, temp_int, temp_int_2)
    ELSE
        ALLOCATE(self % ns_cx(self % max_number_of_nodes))
        ALLOCATE(self % ns_cy(self % max_number_of_nodes))
        ALLOCATE(self % ns_cz(self % max_number_of_nodes))
        ALLOCATE(self % ns_halfside(self % max_number_of_nodes))
        ALLOCATE(self % ns_L2(self % max_number_of_nodes))
        ALLOCATE(self % ns_mass(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_x(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_y(self % max_number_of_nodes))
        ALLOCATE(self % ns_qcm_z(self % max_number_of_nodes))
        ALLOCATE(self % ns_particle(self % max_number_of_nodes))
        ALLOCATE(self % ns_type(self % max_number_of_nodes))
        self % ns_type = 0
        ALLOCATE(self % ns_depth(self % max_number_of_nodes))
        ALLOCATE(self % ns_child(self % max_number_of_nodes, 8))

        ! quadrupole
        ALLOCATE(self % ns_quad(self % max_number_of_nodes, 6))
        self % ns_quad = 0.0_pf
    ENDIF
END SUBROUTINE

FUNCTION node_size_center (x, y, z) RESULT (infos)
! this subroutine gets the node size and center (in space xyz) based in the position
! vectors, giving the minimal square that contains every body.
! @calledby init
    REAL(pf), INTENT(IN) :: x(:), y(:), z(:)
    REAL(pf) :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(pf) :: infos(4)

    xmin = MINVAL(x)
    xmax = MAXVAL(x)
    ymin = MINVAL(y)
    ymax = MAXVAL(y)
    zmin = MINVAL(z)
    zmax = MAXVAL(z)

    infos(1) = 0.5_pf * (xmin + xmax)
    infos(2) = 0.5_pf * (ymin + ymax)
    infos(3) = 0.5_pf * (zmin + zmax)
    infos(4) = MAXVAL((/ xmax - xmin, ymax - ymin, zmax - zmin /))
END FUNCTION

SUBROUTINE add_node (self, cx, cy, cz, side, depth, idx)
! this subroutine adds a node to the tree as a leaf. if its necessary, it reallocate the
! global node state vectors with twice the original size.
! @calledby init, allocate_subnode
    CLASS(OctreeType), INTENT(INOUT) :: self
    REAL(pf), INTENT(IN) :: cx, cy, cz, side
    INTEGER, INTENT(IN) :: depth
    INTEGER, INTENT(INOUT) :: idx

    self % number_of_nodes = self % number_of_nodes + 1
    idx = self % number_of_nodes
    
    IF (idx > self % max_number_of_nodes) THEN
        CALL self % allocate_nodes()
        PRINT *, '[DEBUG] REALLOCATING!'
    END IF

    ! starts without children and being a leaf
    self % ns_child(idx,:) = -1
    self % ns_type(idx) = 1
    self % ns_depth(idx) = depth

    self % ns_cx(idx) = cx
    self % ns_cy(idx) = cy
    self % ns_cz(idx) = cz
    self % ns_halfside(idx) = side / 2.0_pf
    self % ns_L2(idx) = side * side
    self % ns_mass(idx) = 0.0_pf
    self % ns_qcm_x(idx) = 0.0_pf
    self % ns_qcm_y(idx) = 0.0_pf
    self % ns_qcm_z(idx) = 0.0_pf
    self % ns_particle(idx) = -1
    
    ! quadrupole
    self % ns_quad(idx,:) = 0.0_pf
END SUBROUTINE

SUBROUTINE allocate_subnode (self, node_idx, index)
! this subroutine allocates a subnode based in the index of the subnode wrt the node.
! the index is the default: 1-NNO, 2-NNE, 3-NSO, 4-NSE, 5-SNO, 6-SNE, 7-SSO, 8-SSE.
! @calledby add_to_subnode
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: index
    INTEGER :: subnode_idx, d
    REAL(pf) :: h, h_half, cx, cy, cz, cx_sub, cy_sub, cz_sub
    INTEGER :: sign_x, sign_y, sign_z
    
    h = self % ns_halfside(node_idx)
    h_half = h / 2.0_pf
    cx = self % ns_cx(node_idx)
    cy = self % ns_cy(node_idx)
    cz = self % ns_cz(node_idx)
    d = self % ns_depth(node_idx) + 1

    sign_x = MERGE(1,-1,BTEST(index-1,0))
    sign_y = MERGE(-1,1,BTEST(index-1,1))
    sign_z = MERGE(-1,1,BTEST(index-1,2))

    cx_sub = cx + sign_x*h_half
    cy_sub = cy + sign_y*h_half
    cz_sub = cz + sign_z*h_half

    ! create subnode
    CALL self % add_node(cx_sub, cy_sub, cz_sub, h, d, subnode_idx)
    ! allocate it as child of node
    self % ns_child(node_idx, index) = subnode_idx

    IF (self % save_txt .NE. -1) WRITE (self % save_txt, *) d, cx_sub, cy_sub, cz_sub
END SUBROUTINE

SUBROUTINE add_to_subnode (self, node_idx, p)
! given a node/twig and a particle, this adds the particle to the correct leaf based
! in its position.
! @calledby add
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: p
    INTEGER :: subnode
    INTEGER :: subnode_index

    INTEGER :: ix, iy, iz

    ix = MERGE(1,0,self % x(p) > self % ns_cx(node_idx))
    iy = MERGE(1,0,self % y(p) < self % ns_cy(node_idx))
    iz = MERGE(1,0,self % z(p) < self % ns_cz(node_idx))

    subnode_index = 1 + ix + 2*iy + 4*iz
    ! CALL self % index_subnode(node_idx, p, subnode_index)

    subnode = self % ns_child(node_idx,subnode_index)

    ! if the subnode isnt associated
    IF (subnode == -1) THEN
        CALL self % allocate_subnode(node_idx, subnode_index)
    ENDIF

    subnode_index = self % ns_child(node_idx,subnode_index)

    ! now add
    CALL self % add(subnode_index, p)
END SUBROUTINE

SUBROUTINE add (self, node_idx, p)
! given a node and a particle, it adds the particle to the node as a particle if the
! node is empty and as a leaf if the node is already a leaf/particle.
! @calledby init
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node_idx
    INTEGER, INTENT(IN) :: p ! particle index
    INTEGER :: old_p
    REAL(pf) :: pm, px, py, pz, old_mass
    
    ! get particle information
    pm = self % m(p)
    px = self % x(p)
    py = self % y(p)
    pz = self % z(p)

    ! an empty node become a particle
    IF (self % ns_particle(node_idx) == -1 .AND. self % ns_type(node_idx) == 1) THEN
        self % ns_particle(node_idx) = p

        ! add directly the mass and center of mass
        self % ns_mass(node_idx) = pm
        self % ns_qcm_x(node_idx) = px
        self % ns_qcm_y(node_idx) = py
        self % ns_qcm_z(node_idx) = pz

        RETURN
    ENDIF

    ! update the mass and center of mass if the node isnt empty
    old_mass = self % ns_mass(node_idx)
    self % ns_mass(node_idx) = self % ns_mass(node_idx) + pm
    self % ns_qcm_x(node_idx) = (self % ns_qcm_x(node_idx) * old_mass + px * pm) / self % ns_mass(node_idx)
    self % ns_qcm_y(node_idx) = (self % ns_qcm_y(node_idx) * old_mass + py * pm) / self % ns_mass(node_idx)
    self % ns_qcm_z(node_idx) = (self % ns_qcm_z(node_idx) * old_mass + pz * pm) / self % ns_mass(node_idx)

    ! if isnt empty, it become a twig
    IF (self % ns_type(node_idx) == 1) THEN
        ! we cannot go beyond the depth limit
        IF (self % ns_depth(node_idx) >= self % max_depth) THEN
            PRINT *, "BIG PROBLEM !!! MAX DEPTH !!!"
            STOP 0
        ENDIF

        ! in this case, its now a twig
        self % ns_type(node_idx) = 2

        ! add the old particle as a particle per si
        old_p = self % ns_particle(node_idx)
        self % ns_particle(node_idx) = -1
        CALL self % add_to_subnode(node_idx, old_p)
    ENDIF

    ! now add the new particle
    CALL self % add_to_subnode(node_idx, p)
END SUBROUTINE

SUBROUTINE evaluate_quad_aprox (self, par_node_idx)
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, OPTIONAL :: par_node_idx
    INTEGER :: node_idx, child_idx, i
    REAL(pf) :: pm, px, py, pz
    REAL(pf) :: dxi, dyi, dzi ! quadrupole

    node_idx = 1
    IF (PRESENT(par_node_idx)) node_idx = par_node_idx

    IF (node_idx == 1) self % ns_quad = 0.0_pf

    DO i = 1, 8
        child_idx = self % ns_child(node_idx, i)
        
        IF (child_idx == -1) CYCLE

        ! if its a node, evaluate the node
        IF (self % ns_type(child_idx) == 2) THEN
            CALL self % evaluate_quad_aprox(child_idx)
        ENDIF

        ! get information
        pm = self % ns_mass(child_idx)
        px = self % ns_qcm_x(child_idx)
        py = self % ns_qcm_y(child_idx)
        pz = self % ns_qcm_z(child_idx)

        ! update the quadrupole vector
        dxi = px - self % ns_qcm_x(node_idx)
        dyi = py - self % ns_qcm_y(node_idx)
        dzi = pz - self % ns_qcm_z(node_idx)

        self % ns_quad(node_idx, 1) = self % ns_quad(node_idx, 1) + pm * dxi**2    ! mxi2
        self % ns_quad(node_idx, 2) = self % ns_quad(node_idx, 2) + pm * dyi**2    ! myi2
        self % ns_quad(node_idx, 3) = self % ns_quad(node_idx, 3) + pm * dzi**2    ! mzi2
        self % ns_quad(node_idx, 4) = self % ns_quad(node_idx, 4) + pm * dxi * dyi ! mxyi
        self % ns_quad(node_idx, 5) = self % ns_quad(node_idx, 5) + pm * dxi * dzi ! mxzi
        self % ns_quad(node_idx, 6) = self % ns_quad(node_idx, 6) + pm * dyi * dzi ! myzi
    END DO
END SUBROUTINE


SUBROUTINE evaluate_quad (self)
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER :: i, child_idx, p, node_idx
    REAL(pf) :: pm, px, py, pz
    REAL(pf) :: dxi, dyi, dzi ! quadrupole
    INTEGER :: stack(self % number_of_nodes), top

    self % ns_quad = 0.0_pf

    !$OMP PARALLEL SHARED(self) PRIVATE(node_idx, top, stack, child_idx, p, pm, px, py, pz, dxi, dyi, dzi)
    !$OMP DO
    DO node_idx = 1, self % number_of_nodes
        IF (self % ns_type(node_idx) .NE. 2) CYCLE

        top = 0
        DO i = 1, 8
            child_idx = self % ns_child(node_idx, i)

            IF (child_idx .NE. -1) THEN
                top = top + 1
                stack(top) = child_idx
            ENDIF
        END DO
        
        DO WHILE (top > 0)
            p = stack(top)
            top = top - 1

            ! If its a particle, udpate the quadrupole vector
            IF (self % ns_type(p) == 1 .OR. self % ns_depth(p) - self % ns_depth(node_idx) >= 4) THEN
                ! get information
                pm = self % ns_mass(p)
                px = self % ns_qcm_x(p)
                py = self % ns_qcm_y(p)
                pz = self % ns_qcm_z(p)

                ! update the quadrupole vector
                dxi = px - self % ns_qcm_x(node_idx)
                dyi = py - self % ns_qcm_y(node_idx)
                dzi = pz - self % ns_qcm_z(node_idx)

                self % ns_quad(node_idx, 1) = self % ns_quad(node_idx, 1) + pm * dxi**2    ! mxi2
                self % ns_quad(node_idx, 2) = self % ns_quad(node_idx, 2) + pm * dyi**2    ! myi2
                self % ns_quad(node_idx, 3) = self % ns_quad(node_idx, 3) + pm * dzi**2    ! mzi2
                self % ns_quad(node_idx, 4) = self % ns_quad(node_idx, 4) + pm * dxi * dyi ! mxyi
                self % ns_quad(node_idx, 5) = self % ns_quad(node_idx, 5) + pm * dxi * dzi ! mxzi
                self % ns_quad(node_idx, 6) = self % ns_quad(node_idx, 6) + pm * dyi * dzi ! myzi
            
            ! If its a twig, update the twig
            ELSE IF (self % ns_type(p) == 2) THEN
                ! then add the children to the stack
                DO i = 1, 8
                    child_idx = self % ns_child(p, i)
                    IF (child_idx .NE. -1) THEN
                        top = top + 1
                        stack(top) = child_idx
                    ENDIF
                END DO
            ENDIF
        END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
END SUBROUTINE

FUNCTION evaluate_forces_over_p (self, p, par_theta2, par_G, par_eps2) RESULT (forces)
! given a particle and the parameters, this evaluates the forces over the particle
! using the Barnes-Hut criterion. it uses the depth first traversal (DFS) algorithm 
! to evaluate the forces in the tree.
    CLASS(OctreeType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: p ! particle index
    REAL(pf), INTENT(IN), OPTIONAL :: par_theta2, par_G, par_eps2 ! parameters
    REAL(pf) :: eps2, theta2, G

    REAL(pf) :: pm, px, py, pz
    REAL(pf) :: forces(3)
    REAL(pf) :: dx, dy, dz, dist2, L2, f

    INTEGER :: stack(self % number_of_nodes)
    INTEGER :: top, node_idx, i, child_idx

    REAL(pf) :: rinv, invR, invR3, invR5
    REAL(pf) :: Mx, My, Mz
    REAL(pf) :: Qx, Qy, Qz
    REAL(pf) :: mxi2, myi2, mzi2, mxyi, mxzi, myzi

    REAL(pf) :: dx2, dy2, dz2, dxyz, invR2

    ! default values
    eps2 = 0.0_pf
    theta2 = 0.0_pf
    G = 1.0_pf

    ! replace if present
    IF (PRESENT(par_eps2))   eps2 = par_eps2
    IF (PRESENT(par_theta2)) theta2 = par_theta2
    IF (PRESENT(par_G))      G = par_G

    ! particle cache info
    pm = self % m(p)
    px = self % x(p)
    py = self % y(p)
    pz = self % z(p)

    ! initialize
    forces = 0.0_pf
    top = 1
    stack(top) = 1

    ! dfs iterative
    DO WHILE (top > 0)

        node_idx = stack(top)
        top = top - 1

        ! empty node
        IF (self % ns_mass(node_idx) == 0.0_pf) CYCLE

        ! self interaction
        IF (self % ns_type(node_idx) == 1) THEN
            IF (self % ns_particle(node_idx) == p) CYCLE
        ENDIF

        ! geometry
        dx = self % ns_qcm_x(node_idx) - px
        dy = self % ns_qcm_y(node_idx) - py
        dz = self % ns_qcm_z(node_idx) - pz
        dist2 = dx*dx + dy*dy + dz*dz
        L2 = self % ns_L2(node_idx)

        ! leaf or bh criterion
        IF (self % ns_type(node_idx) == 1 .OR. L2 < theta2 * dist2) THEN
            IF (self % ns_type(node_idx) == 1) THEN
                rinv = 1.0_pf / SQRT(dist2 + eps2)
                rinv = rinv * rinv * rinv
                f = G * pm * self % ns_mass(node_idx) * rinv
                forces(1) = forces(1) + f * dx
                forces(2) = forces(2) + f * dy
                forces(3) = forces(3) + f * dz
            ELSE
                invR = 1.0_pf / SQRT(dist2 + eps2)
                invR2 = 1.0_pf / (dist2 + eps2)
                invR3 = invR * invR2
                invR5 = invR3 * invR * invR

                Mx = dx * invR3 * self % ns_mass(node_idx)
                My = dy * invR3 * self % ns_mass(node_idx)
                Mz = dz * invR3 * self % ns_mass(node_idx)
                
                mxi2 = self % ns_quad(node_idx, 1)
                myi2 = self % ns_quad(node_idx, 2)
                mzi2 = self % ns_quad(node_idx, 3)
                mxyi = self % ns_quad(node_idx, 4)
                mxzi = self % ns_quad(node_idx, 5)
                myzi = self % ns_quad(node_idx, 6)

                dx2  = 15.0_pf * dx * dx * invR2 - 3.0_pf
                dy2  = 15.0_pf * dy * dy * invR2 - 3.0_pf
                dz2  = 15.0_pf * dz * dz * invR2 - 3.0_pf
                dxyz = 15.0_pf * dx * dy * dz * invR2
                
                Qx = (mxi2 * (dx2 - 6.0_pf) + myi2 * dy2 + mzi2 * dz2) * 0.5_pf * dx
                Qy = (myi2 * (dy2 - 6.0_pf) + mzi2 * dz2 + mxi2 * dx2) * 0.5_pf * dy
                Qz = (mzi2 * (dz2 - 6.0_pf) + mxi2 * dx2 + myi2 * dy2) * 0.5_pf * dz

                Qx = Qx + mxyi * dy * dx2 + mxzi * dz * dx2 + myzi * dxyz
                Qy = Qy + myzi * dz * dy2 + mxyi * dx * dy2 + mxzi * dxyz
                Qz = Qz + mxzi * dx * dz2 + myzi * dy * dz2 + mxyi * dxyz

                forces(1) = forces(1) + G * pm * (Mx + invR5 * Qx)
                forces(2) = forces(2) + G * pm * (My + invR5 * Qy)
                forces(3) = forces(3) + G * pm * (Mz + invR5 * Qz)
            ENDIF
        ! if isnt leaf nor bh is valid, so go to children
        ELSE
            ! push children in the vec
            DO i = 1, 8
                child_idx = self % ns_child(node_idx, i)
                IF (child_idx .NE. -1) THEN
                    top = top + 1
                    stack(top) = child_idx
                ENDIF
            END DO
        ENDIF
    END DO
END FUNCTION

END MODULE 