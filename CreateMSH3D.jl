function CreateMSH3D(name)

    gmsh.initialize()

    if name == "Cube"

        gmsh.model.add("Cube")
        lc = 0.1
        p1 = gmsh.model.geo.addPoint(-1, -1, -1, lc)
        p2 = gmsh.model.geo.addPoint(1, -1, -1, lc)
        p3 = gmsh.model.geo.addPoint(1, 1, -1, lc)
        p4 = gmsh.model.geo.addPoint(-1, 1, -1, lc)
        p5 = gmsh.model.geo.addPoint(-1, -1, 1, lc)
        p6 = gmsh.model.geo.addPoint(1, -1, 1, lc)
        p7 = gmsh.model.geo.addPoint(1, 1, 1, lc)
        p8 = gmsh.model.geo.addPoint(-1, 1, 1, lc)

        l1 = gmsh.model.geo.addLine(p1, p2)
        l2 = gmsh.model.geo.addLine(p2, p3)
        l3 = gmsh.model.geo.addLine(p3, p4)
        l4 = gmsh.model.geo.addLine(p4, p1)
        l5 = gmsh.model.geo.addLine(p5, p6)
        l6 = gmsh.model.geo.addLine(p6, p7)
        l7 = gmsh.model.geo.addLine(p7, p8)
        l8 = gmsh.model.geo.addLine(p8, p5)
        l9 = gmsh.model.geo.addLine(p1, p5)
        l10 = gmsh.model.geo.addLine(p2, p6)
        l11 = gmsh.model.geo.addLine(p3, p7)
        l12 = gmsh.model.geo.addLine(p4, p8)

        cl1 = gmsh.model.geo.addCurveLoop([-l1, -l4, -l3, -l2])
        cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
        cl3 = gmsh.model.geo.addCurveLoop([l1,l10,-l5,-l9])
        cl4 = gmsh.model.geo.addCurveLoop([l2,l11,-l6,-l10])
        cl5 = gmsh.model.geo.addCurveLoop([l3,l12,-l7,-l11])
        cl6 = gmsh.model.geo.addCurveLoop([l4,l9,-l8,-l12])

        ps1 = gmsh.model.geo.addPlaneSurface([cl1])
        ps2 = gmsh.model.geo.addPlaneSurface([cl2])
        ps3 = gmsh.model.geo.addPlaneSurface([cl3])
        ps4 = gmsh.model.geo.addPlaneSurface([cl4])
        ps5 = gmsh.model.geo.addPlaneSurface([cl5])
        ps6 = gmsh.model.geo.addPlaneSurface([cl6])

        s1 = gmsh.model.geo.addSurfaceLoop([ps1, ps2, ps3, ps4, ps5, ps6])

        v1 = gmsh.model.geo.addVolume([s1])

        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)

        _, POS, _ = gmsh.model.mesh.getNodes()
        POS = reshape(POS, (3,:))'
        N_nodes = size(POS,1)

        _, _, TETRAHEDRA = gmsh.model.mesh.getElements(3)
        TETRAHEDRA = reshape(TETRAHEDRA[1], (4,:))'
        N_tetra = size(TETRAHEDRA,1)

        directions = ["D", "U", "F", "R", "B", "L"]
        TRIANGLES = Dict()
        for (ps_i, dir) ∈ enumerate(directions)
            _, _, TRIANGLES[dir] = gmsh.model.mesh.getElements(2, ps_i)
            TRIANGLES[dir]= reshape(TRIANGLES[dir][1], (3,:))'
        end

        gmsh.write("xx_msh/Cube.msh")

    end

    gmsh.finalize()

    return N_nodes, N_tetra, POS, TETRAHEDRA, TRIANGLES

end
