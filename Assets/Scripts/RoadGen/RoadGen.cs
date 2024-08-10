using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using UnityEngine.Rendering;

public static class RoadGen
{    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="spline"></param>
    /// <param name="partMesh"></param>
    /// <param name="roadWidth">Width of the road. If -1, it should be calculated automatically</param>
    /// <param name="roadPartLength">Length of segment. Note that this is only a rough value, as it will be adapted to be a divisor of the total curve length, in order not to end the loop somewhere in the segmentMesh</param>
    /// <param name="roadPartHeight">Height of one part of the road</param>
    public static void GenerateRoad_WholeParts(Transform roadObject, RoadSpline spline, float widthMultiplier, float heightMultiplier, float crossingShape, Dictionary<TopologyType, MeshConfig> meshConfigs)
    {
        var topologies = TopologyHandler.topologies;
        TopologyHandler.LogTopologies(topologies);
        if (topologies == null)
        {
            Debug.LogWarning("ERROR: Must generate topologies first.");
            return;
        }

        //https://docs.unity3d.com/ScriptReference/Object.DestroyImmediate.html
        var children = roadObject.Cast<Transform>().ToList();
        foreach (var child in children)
        {
            UnityEngine.Object.DestroyImmediate(child.gameObject);
        }


        bool hasWrappingTopology = topologies[0].type == topologies[^1].type && topologies.Count > 1; // Practically, this should always be true when having more than one topology

        #region Everything but Crossings and Wrapping Topologies

        {
            float lastPercEnd = 0;
            for (int i = 0; i < topologies.Count; i++)
            {
                TopologyRange topologyRange = topologies[i];
                TopologyType topologyType = topologyRange.type;
        
                float percEnd = topologyRange.percEnd;
                bool wrappingTopology = (i == 0 || i == topologies.Count - 1) && hasWrappingTopology;
        
                if (topologyType != TopologyType.Crossing && !wrappingTopology)
                {
        
                    Mesh partMesh = meshConfigs[topologyType].partMesh;
                    float partLength = meshConfigs[topologyType].partLength;
        
                    Mesh roadMesh = CreateMeshForTopology(
                        spline: spline,
                        percStart: lastPercEnd,
                        percEnd: percEnd,
                        partLength: partLength,
                        widthMultiplier: widthMultiplier,
                        heightMultiplier: heightMultiplier,
                        partMesh: partMesh,
                        out _
                    );
        
                    List<Material> materials = meshConfigs[topologyType].materials;
                    string name = topologyRange.type.ToString();// + "_" + Time.time;
        
                    MeshAsGameObject(mesh: roadMesh, materials: materials, name: name).transform.SetParent(roadObject);
                }
                lastPercEnd = percEnd;
            }
        }

        #endregion

        #region Wrapping Topology

        if (hasWrappingTopology)
        {
            TopologyType topologyType = topologies[0].type;
            if (topologyType == TopologyType.Crossing)
            {
                Debug.LogWarning("ERROR: Wrapping Crossing not supported yet -- sorry.");
            }
            else
            {
                float percStart = topologies[^2].percEnd;
                float percEnd = topologies[0].percEnd + 1f;

                Mesh partMesh = meshConfigs[topologyType].partMesh;
                float partLength = meshConfigs[topologyType].partLength;

                Mesh roadMesh = CreateMeshForTopology(
                    spline: spline,
                    percStart: percStart,
                    percEnd: percEnd,
                    partLength: partLength,
                    widthMultiplier: widthMultiplier,
                    heightMultiplier: heightMultiplier,
                    partMesh: partMesh,
                    out _
                );

                List<Material> materials = meshConfigs[topologyType].materials;
                string name = topologies[0].type.ToString() + "_" + Time.time;

                MeshAsGameObject(mesh: roadMesh, materials: materials, name: name).transform.SetParent(roadObject);
            }
        }

        #endregion

        #region Crossings

        List<List<float>> crossings = TopologyHandler.fourPercCrossings;
        foreach (var c in crossings) // Crossings
        {
            Vector3 roughMidpoint = Vector3.zero;
            for (int q = 0; q < c.Count; q++) // Quarters
            {
                float a_perc = c[q];
                float b_perc = c[(q + 1) % c.Count];

                float a_dst = a_perc * spline.TotalLength;
                float b_dst = b_perc * spline.TotalLength;

                var (a_anchorPoint, a_segIndex, a_t) = spline.CalculatePointAtDist(a_dst);
                var (b_anchorPoint, b_segIndex, b_t) = spline.CalculatePointAtDist(b_dst);

                Vector2 a_tangent = spline.UnitTangent(a_segIndex, a_t);
                Vector2 b_tangent = spline.UnitTangent(b_segIndex, b_t);

                // Note that the Intersection of a_tangent and b_tangent is NOT the global intersection point of the curve, as the crossing sections are not linear
                var (a_ctrlMaxLength, b_ctrlMaxLength) = MyMath.Intersection(a_anchorPoint, a_tangent, b_anchorPoint, b_tangent);

                Vector2 tangentIntersectionPoint = a_anchorPoint + a_tangent * a_ctrlMaxLength;
                float a_ctrlLength = a_ctrlMaxLength * crossingShape;
                float b_ctrlLength = b_ctrlMaxLength * crossingShape;

                Vector2 a_ctrlPoint = a_anchorPoint + a_tangent * a_ctrlLength;
                Vector2 b_ctrlPoint = b_anchorPoint + b_tangent * b_ctrlLength;

                RoadSpline quarterSpline = new(
                    new List<Vector2>() { a_anchorPoint, a_ctrlPoint, b_ctrlPoint, b_anchorPoint },
                    true,
                    false
                );

                TopologyType topologyType = TopologyType.Crossing;

                float partLength = meshConfigs[topologyType].partLength;
                Mesh partMesh = meshConfigs[topologyType].partMesh;

                Mesh quarterMesh = CreateMeshForTopology(
                    spline: quarterSpline,
                    percStart: 0,
                    percEnd: 1,
                    partLength: partLength,
                    widthMultiplier: widthMultiplier,
                    heightMultiplier: heightMultiplier,
                    partMesh: partMesh,
                    out List<List<int>> indexesOfMidVerts
                );

                List<Vector3> vertices = quarterMesh.vertices.ToList();
                if (q == 0)
                {
                    int firstIndex = indexesOfMidVerts.SelectMany(list => list).ToList()[0];
                    roughMidpoint = new Vector3(tangentIntersectionPoint.x, vertices[firstIndex].y, tangentIntersectionPoint.y);
                }
                vertices.Add(roughMidpoint);
                int intersectionVertIndex = vertices.Count - 1;
                quarterMesh.vertices = vertices.ToArray(); // Must be set before we do SetTriangles()

                for (int sm = 0; sm < quarterMesh.subMeshCount; sm++)
                {
                    var smTriangles = quarterMesh.GetTriangles(sm).ToList();
                    for (int i = 0; i < indexesOfMidVerts[sm].Count; i += 2)
                    {
                        smTriangles.AddRange(new List<int>() {
                            indexesOfMidVerts[sm][i],
                            indexesOfMidVerts[sm][i + 1],
                            intersectionVertIndex
                        });
                    }
                    quarterMesh.SetTriangles(smTriangles, sm);
                }
                quarterMesh.RecalculateBounds();
                quarterMesh.RecalculateNormals();

                List<Material> materials = meshConfigs[topologyType].materials;
                string name = "CrossingQuarter" + q + "_" + Time.time;

                MeshAsGameObject(mesh: quarterMesh, materials: materials, name: name).transform.SetParent(roadObject);
            }
        }

        #endregion



        #region DEBUG

        //bool firstSegmentDone = false;
        //float ref_totalDstAlongSpline = 0;

        //float DEBUG_LENGTH_MULTIPLIER = 1f;//1.57f;


        //Debug.Log("ref_totalDstAlongSpline at End: " + ref_totalDstAlongSpline);
        //Debug.Log("but total length is: " + totalLength);
        //Debug.Log("Length by parts is: " + (numParts * meshLength * meshLengthToRoadLength));
        //Debug.Log("numParts: " + numParts);
        //Debug.Log("partLength: " + roadPartLength);
        //Debug.Log("meshLength * meshLengthToRoadLength: " + meshLength * meshLengthToRoadLength);

        {
            //RoadSpline testSpline = new(new List<Vector2>() { new(0f, 0f), new(0.25f, 0f), new(0.75f, 0f), new(1f, 0f) }, true);
            //int segIndex = 0;
            //int refT = 0;
            //float dstAfterRef = 0.5f;
            //float t = testSpline.CalculatePointAfterRef(dstAfterRef, segIndex, refT, 10f).t;
            //Debug.Log(testSpline);
            //Debug.Log("(wrong) t for Dst of " + dstAfterRef + " is " + t);
            //Debug.Log("(wrong) dst on refPoint is: " + testSpline.GetDstAlongSpline(segIndex, refT));
            //Debug.Log("(wrong) dst for newPoint is: " + testSpline.GetDstAlongSpline(segIndex, t));

            //Debug.Log(spline.CalculatePointAfterRef(spline.TotalLength, 0, 0));
            //Debug.Log(spline.CalculatePointAtDist(spline.TotalLength, 10f));
            //Debug.Log(spline.FirstPoint);
            //Debug.Log(spline.GetLastPoint().point);
        }

        #endregion
    }

    static Mesh CreateMeshForTopology(RoadSpline spline, float percStart, float percEnd, float partLength, float widthMultiplier, float heightMultiplier, Mesh partMesh, out List<List<int>> indexesOfMidVerts)
    {
        /* Naming Convention:
         * mesh...: Refers to the mesh of the road part that was given as an argument
         * ref...: Refers to the ref vert that is used for calculating the next points from a known point rather than from the start
         * spline...: Refers to the Vertex etc. on the spline
         * road...: Refers to (part of) the road mesh that is being generated
         */

        float splineLength = spline.TotalLength;
        float topologyLength = splineLength * (percEnd - percStart);
        int numParts = (int)(topologyLength / partLength);
        numParts = Mathf.Max(numParts, 1); // Two reasons: (1) If topologyLength is roughly partLength, numParts can be zero. (2) If a segment ends right before t = 1, then numParts could also be zero
        partLength = topologyLength / numParts;

        var meshVertices = partMesh.vertices.ToList();//.OrderBy(x => x.z).ToList(); // Ordering by z (length) is helpful for calculating the positions on the spline.
        float meshMaxZ = meshVertices.ToList().OrderBy(v => v.z).Last().z;
        float meshMinZ = meshVertices.ToList().OrderBy(v => v.z).First().z;
        float meshWidth = partMesh.bounds.size.x;
        float meshHeight = partMesh.bounds.size.y;
        float meshLength = partMesh.bounds.size.z;
        float meshLengthToRoadLength = partLength / meshLength;

        var vertices = new List<Vector3>();
        var triangles = new List<List<int>>(); // List of List of triangles of every submesh

        float PercInMesh(Vector3 meshVert) => (meshVert.z - meshMinZ) / meshLength;

        indexesOfMidVerts = new();

        List<List<int>> trianglesOfSubMeshes = new(); // Here for garbage collection
        for (int sm = 0; sm < partMesh.subMeshCount; sm++)
        {
            trianglesOfSubMeshes.Add(partMesh.GetTriangles(sm).ToList());
        }

        // Loop through all parts
        for (int p = 0; p < numParts; p++)
        {
            // Loop through subMeshes and add Triangles
            for (int sm = 0; sm < partMesh.subMeshCount; sm++)
            {
                int[] smTriangles = new List<int>(trianglesOfSubMeshes[sm]).ToArray();
                for (int tri = 0; tri < smTriangles.Length; tri++)
                    smTriangles[tri] += vertices.Count;
                if (sm >= triangles.Count) triangles.Add(new());
                triangles[sm].AddRange(smTriangles);
            }

            // Loop through subMeshes and add vertices
            Vector3 spline_point = Vector3.zero; int spline_segIndex = 0; float spline_t = 0; // Needed here to override ref
            for (int sm = 0; sm < partMesh.subMeshCount; sm++)
            {
                SubMeshDescriptor subMesh = partMesh.GetSubMesh(sm);

                for (int v = 0; v < subMesh.vertexCount; v++)
                {
                    // Get spacing along spline in road coordinates to ref point
                    Vector3 meshVert = meshVertices[v + subMesh.firstVertex];
                    float meshVert_z = meshVert.z;

                    // Get the 3D position of the vert on the spline
                    (spline_point, spline_segIndex, spline_t) = spline.CalculatePointAtDist((percStart * splineLength) + (p * partLength) + (PercInMesh(meshVert) * partLength));
                    Vector3 splineVert_3D = new(spline_point.x, 0, spline_point.y);                 // 2D to 3D... only for now, until Spline operates in 3D


                    // Get the normal of the vert on the spline
                    Vector3 splineVert_tangent = spline.UnitTangent(spline_segIndex, spline_t);
                    splineVert_tangent = new(splineVert_tangent.x, 0, splineVert_tangent.y);    // 2D to 3D... only for now, until Spline operates in 3D
                    Vector3 splineVert_normal = Vector3.Cross(splineVert_tangent, Vector3.up);

                    // Get the x offset of the vert
                    float roadVert_offsetWidth = meshVert.x * widthMultiplier;

                    // Get the y offset of the vert
                    float roadVert_offsetHeight = meshVert.y * heightMultiplier;

                    if (sm >= indexesOfMidVerts.Count) indexesOfMidVerts.Add(new());
                    if (Mathf.Abs(meshVert.x) < 0.0001f)
                        indexesOfMidVerts[sm].Add(vertices.Count);

                    Vector3 roadVert = splineVert_3D + splineVert_normal * roadVert_offsetWidth + Vector3.up * roadVert_offsetHeight;
                    vertices.Add(roadVert);
                }
            }

            //ref_totalDstAlongSpline += roadPartLength * DEBUG_LENGTH_MULTIPLIER;
        }

        Mesh mesh = new()
        {
            subMeshCount = partMesh.subMeshCount
        };

        mesh.SetVertices(vertices);
        for (int i = 0; i < partMesh.subMeshCount; i++)
            mesh.SetTriangles(triangles[i], i);
        mesh.RecalculateBounds();
        mesh.RecalculateNormals();

        return mesh;
    }

    static Mesh CreateMeshForTopology_UsingRef(RoadSpline spline, float percStart, float percEnd, float partLength, float widthMultiplier, float heightMultiplier, Mesh partMesh, out List<List<int>> indexesOfMidVerts)
    {
        /* Naming Convention:
         * mesh...: Refers to the mesh of the road part that was given as an argument
         * ref...: Refers to the ref vert that is used for calculating the next points from a known point rather than from the start
         * spline...: Refers to the Vertex etc. on the spline
         * road...: Refers to (part of) the road mesh that is being generated
         */

        float splineLength = spline.TotalLength;
        float topologyLength = splineLength * (percEnd - percStart);
        int numParts = (int)(topologyLength / partLength);
        numParts = Mathf.Max(numParts, 1); // Two reasons: (1) If topologyLength is roughly partLength, numParts can be zero. (2) If a segment ends right before t = 1, then numParts could also be zero
        partLength = topologyLength / numParts;

        var meshVertices = partMesh.vertices.ToList();//.OrderBy(x => x.z).ToList(); // Ordering by z (length) is helpful for calculating the positions on the spline.
        float meshMaxZ = meshVertices.ToList().OrderBy(v => v.z).Last().z;
        float meshMinZ = meshVertices.ToList().OrderBy(v => v.z).First().z;
        float meshWidth = partMesh.bounds.size.x;
        float meshHeight = partMesh.bounds.size.y;
        float meshLength = partMesh.bounds.size.z;
        float meshLengthToRoadLength = partLength / meshLength;

        var vertices = new List<Vector3>();
        var triangles = new List<List<int>>(); // List of List of triangles of every submesh

        (Vector3 ref_spline_point, int ref_spline_segIndex, float ref_spline_t) = spline.CalculatePointAtDist(percStart * spline.TotalLength);
        float ref_meshVert_z = meshMinZ;

        indexesOfMidVerts = new();

        // Loop through all parts
        for (int p = 0; p < numParts; p++)
        {
            // Loop through subMeshes and add Triangles
            for (int sm = 0; sm < partMesh.subMeshCount; sm++)
            {
                var smTriangles = partMesh.GetTriangles(sm);
                for (int tri = 0; tri < smTriangles.Length; tri++)
                    smTriangles[tri] += vertices.Count;
                if (sm >= triangles.Count) triangles.Add(new());
                triangles[sm].AddRange(smTriangles);
            }

            var next_ref_spline_point = ref_spline_point;
            var next_ref_spline_segIndex = ref_spline_segIndex;
            var next_ref_spline_t = ref_spline_t;

            // Loop through subMeshes and add vertices
            Vector3 spline_point = Vector3.zero; int spline_segIndex = 0; float spline_t = 0; // Needed here to override ref
            for (int sm = 0; sm < partMesh.subMeshCount; sm++)
            {
                SubMeshDescriptor subMesh = partMesh.GetSubMesh(sm);

                for (int v = 0; v < subMesh.vertexCount; v++)
                {
                    // Get spacing along spline in road coordinates to ref point
                    Vector3 meshVert = meshVertices[v + subMesh.firstVertex];
                    float meshVert_z = meshVert.z;
                    float mesh_spacingToRef = meshVert_z - ref_meshVert_z;
                    float spline_spacingToRef = mesh_spacingToRef * meshLengthToRoadLength;

                    // Get the 3D position of the vert on the spline
                    if (p == numParts - 1 && Mathf.Abs(meshVert_z - meshMaxZ) < 0.0001f) // We're in the very last edgeLoop
                    {
                        (spline_point, spline_segIndex, spline_t) = spline.CalculatePointAtDist(percEnd * spline.TotalLength);
                    }
                    else if (spline_spacingToRef == 0) // Same z as ref_point?
                        (spline_point, spline_segIndex, spline_t) = (ref_spline_point, ref_spline_segIndex, ref_spline_t);
                    else // We need to calculate (usual case)
                        (spline_point, spline_segIndex, spline_t) = spline.CalculatePointAfterRef(spline_spacingToRef, ref_spline_segIndex, ref_spline_t);
                    Vector3 splineVert_3D = new(spline_point.x, 0, spline_point.y);                 // 2D to 3D... only for now, until Spline operates in 3D

                    //Debug.Log($"({spline_t}, {spline_segIndex})");

                    // prepare next ref
                    if (meshVert_z == meshMaxZ)
                    {
                        next_ref_spline_point = spline_point;
                        next_ref_spline_segIndex = spline_segIndex;
                        next_ref_spline_t = spline_t;
                    }

                    // Get the normal of the vert on the spline
                    Vector3 splineVert_tangent = spline.UnitTangent(spline_segIndex, spline_t);
                    splineVert_tangent = new(splineVert_tangent.x, 0, splineVert_tangent.y);    // 2D to 3D... only for now, until Spline operates in 3D
                    Vector3 splineVert_normal = Vector3.Cross(splineVert_tangent, Vector3.up);

                    // Get the x offset of the vert
                    float roadVert_offsetWidth = meshVert.x * widthMultiplier;

                    // Get the y offset of the vert
                    float roadVert_offsetHeight = meshVert.y * heightMultiplier;

                    if (sm >= indexesOfMidVerts.Count) indexesOfMidVerts.Add(new());
                    if (Mathf.Abs(meshVert.x) < 0.0001f)
                        indexesOfMidVerts[sm].Add(vertices.Count);

                    Vector3 roadVert = splineVert_3D + splineVert_normal * roadVert_offsetWidth + Vector3.up * roadVert_offsetHeight;
                    vertices.Add(roadVert);
                }
            }


            ref_spline_point = next_ref_spline_point;
            ref_spline_segIndex = next_ref_spline_segIndex;
            ref_spline_t = next_ref_spline_t;

            //ref_totalDstAlongSpline += roadPartLength * DEBUG_LENGTH_MULTIPLIER;
        }

        Mesh mesh = new()
        {
            vertices = vertices.ToArray(),
            subMeshCount = partMesh.subMeshCount
        };

        for (int i = 0; i < partMesh.subMeshCount; i++)
            mesh.SetTriangles(triangles[i], i);
        mesh.RecalculateBounds();
        mesh.RecalculateNormals();

        return mesh;
    }

    static GameObject MeshAsGameObject(Mesh mesh, List<Material> materials, string name)
    {
        GameObject obj = new()
        {
            name = name
        };

        obj.AddComponent<MeshRenderer>();
        obj.AddComponent<MeshFilter>();
        obj.AddComponent<MeshCollider>();

        obj.GetComponent<MeshRenderer>().materials = materials.ToArray();
        obj.GetComponent<MeshFilter>().sharedMesh = mesh;
        obj.GetComponent<MeshCollider>().sharedMesh = mesh;

        return obj;
    }

    /*

    public static void GenerateRoad_2DSliceExtrude(GameObject roadObject, RoadSpline spline, Mesh partMesh, List<Material> materials, float roadWidth, float roadPartLength, float roadPartHeight, AnimationCurve bridgeStartEndCurve, float bridgeHeight)
    {
        #region Preparation

        List<Vector2> slicePoints = new();

        foreach (var v in partMesh.vertices)
        {
            if (v.z > 0)
            {
                slicePoints.Add(new Vector2(v.x, v.y));
            }
        }

        int edgeLoops = Mathf.RoundToInt(spline.TotalLength / roadPartLength);
        roadPartLength = spline.TotalLength / edgeLoops;
        int vertCount = slicePoints.Count * edgeLoops;
        int numTris = 2 * slicePoints.Count * edgeLoops;
        int numTrisIncides = numTris * 3;

        int[] trisIndices = new int[numTrisIncides];
        Vector3[] vertices = new Vector3[vertCount];
        Vector3[] normals = new Vector3[vertCount];
        Vector2[] uvs = new Vector2[vertCount];

        float meshWidthToRoadWidth = roadWidth / partMesh.bounds.size.x;
        float meshHeightToRoadHeight = roadPartHeight / partMesh.bounds.size.x;

        #endregion

        #region Mesh Generation

        int ref_spline_segIndex = 0;
        float ref_spline_t = 0;
        Vector2 ref_lastPoint = spline.FirstPoint;

        for (int loop = 0; loop < edgeLoops; loop++)
        {
            int vertIndexOffset = loop * slicePoints.Count;
            (Vector2 pointOnSpline, int spline_segIndex, float spline_t) = spline.CalculatePoint(roadPartLength, ref_lastPoint, ref_spline_segIndex, ref_spline_t);
            for (int v = 0; v < slicePoints.Count; v++)
            {
                int id = vertIndexOffset + v;

                Vector3 normal = xyToX0Z(spline.Normal(spline_segIndex, spline_t));

                Vector3 vert = xyToX0Z(pointOnSpline);
                vert += normal * slicePoints[v].x * meshWidthToRoadWidth;
                vert += Vector3.up * slicePoints[v].y * meshHeightToRoadHeight;
                vertices[id] = vert;
                normals[id] = normal;
                uvs[id] = Vector2.zero; // ToDo
            }

            ref_spline_segIndex = spline_segIndex;
            ref_spline_t = spline_t;
            ref_lastPoint = pointOnSpline;
        }

        int ti = 0;
        for (int loop = 1; loop < edgeLoops; loop++)
        {
            int offset_prev = (loop - 1) * slicePoints.Count;
            int offset_curr = loop * slicePoints.Count;
            for (int v = 1; v < slicePoints.Count; v++)
            {
                int a = offset_prev + v - 1;
                int b = offset_prev + v;
                int c = offset_curr + v - 1;
                int d = offset_curr + v;

                trisIndices[ti] = a; ti++;
                trisIndices[ti] = b; ti++;
                trisIndices[ti] = c; ti++;
                trisIndices[ti] = c; ti++;
                trisIndices[ti] = a; ti++;
                trisIndices[ti] = d; ti++;
            }
        }

        Vector3 xyToX0Z(Vector2 xy)
        {
            return new Vector3(xy.x, 0, xy.y);
        }

        #endregion

        #region Apply

        if (roadObject.GetComponent<MeshRenderer>() == null)
            roadObject.AddComponent<MeshRenderer>();
        if (roadObject.GetComponent<MeshFilter>() == null)
            roadObject.AddComponent<MeshFilter>();

        Mesh mesh = roadObject.GetComponent<MeshFilter>().sharedMesh;
        mesh.Clear();
        mesh.vertices = vertices;
        mesh.triangles = trisIndices;
        mesh.normals = normals;
        mesh.uv = uvs;
        mesh.RecalculateBounds();
        mesh.RecalculateNormals();

        roadObject.GetComponent<MeshRenderer>().materials = materials.ToArray();

        #endregion
    }

    public static void GenerateRoad_NowForReal(GameObject roadObject, RoadSpline spline, Mesh roadRef, List<Material> materials, float roadPartWidth, float roadPartLength, float roadPartHeight, AnimationCurve bridgeStartEndCurve, float bridgeHeight)
    {
        int partCount = Mathf.RoundToInt(spline.TotalLength / roadPartLength);
        roadPartLength = spline.TotalLength / partCount;

        List<int> roadRefTriangles = roadRef.triangles.ToList();
        List<Vector2> roadRefUVs = roadRef.uv.ToList();
        List<Vector3> roadRefVerts = roadRef.vertices.ToList();
        float roadRefLength = roadRef.bounds.size.z;
        float roadRefWidth = roadRef.bounds.size.x;
        float roadRefHeight = roadRef.bounds.size.y;
        float roadRefMidY = roadRef.bounds.center.y;

        float esp_spacing = 0.1f;
        spline.CalculateEvenlySpacedPoints(esp_spacing, out List<int> indicesOfSegments, out List<float> tsOfESP);

        MeshWrapper mw = new();

        float dstToLastSeg = 0;

        for (int seg = 0; seg < spline.NumSegments; seg++)
        {
            float segLength = spline.SegmentLength(seg);
            int partsInThisSegment = (int)((segLength - dstToLastSeg) / roadPartLength + 0.001f); // + 0.001f: Making sure we generate the last part of the last segment

            MeshWrapper segMW = new();

            for (int part = 0; part < partsInThisSegment; part++)
            {
                MeshWrapper partMW = AdaptedMesh(dstToLastSeg + roadPartLength * part,
                    roadRefVerts.Count * part + mw.vertsPos.Count,
                    segLength);
                segMW.Add(partMW);
            }

            List<float> distancedTs = DistancedTs(seg, segMW.vertsPercLength);

            List<Vector2> splinePoints2D = spline.CalculatePointsOnCurve(seg, distancedTs.ToArray()).ToList();

            List<Vector2> forwards2D = spline.CalculateForwards(seg, splinePoints2D.ToArray(), distancedTs.ToArray()).ToList();

            float spline_t = 0.5f;
            for (int i = 0; i < splinePoints2D.Count; i++)
            {
                //empty = Instantiate(emptyPref, points[i], Quaternion.identity).transform;
                Vector2 forward2D = forwards2D[i];
                //print(points[i]);
                Vector2 normal2D = new(-forward2D.y, forward2D.x);

                Vector3 currPos = segMW.vertsPos[i];

                Vector3 newPos = new(splinePoints2D[i].x, 0, splinePoints2D[i].y);
                newPos += new Vector3(normal2D.x, 0, normal2D.y) * (segMW.vertsPercWidth[i] * roadPartWidth / 2f);
                newPos += roadPartHeight * segMW.vertsPercHeight[i] * Vector3.up;

                SegmentType segmentType = spline.GetSegmentType(seg);
                if (segmentType == SegmentType.BridgeStart && currPos.y > roadRefMidY)
                    newPos += Vector3.up * bridgeStartEndCurve.Evaluate(spline_t) * bridgeHeight;
                else if (segmentType == SegmentType.Bridge)
                    newPos += Vector3.up * bridgeHeight;
                else if (segmentType == SegmentType.BridgeEnd && currPos.y > roadRefMidY)
                    newPos += Vector3.up * bridgeStartEndCurve.Evaluate(1 - spline_t) * bridgeHeight;

                segMW.vertsPos[i] = newPos;
            }


            dstToLastSeg = dstToLastSeg + roadPartLength * partsInThisSegment - segLength;

            mw.Add(segMW);
        }

        MeshWrapper AdaptedMesh(float startValue, int trisAddition, float segLength)
        {
            MeshWrapper newMesh = new();

            newMesh.tris.AddRange(roadRefTriangles);
            newMesh.uvs = roadRefUVs;

            for (int i = 0; i < roadRefVerts.Count; i++)
            {
                float vertsPercLength = (startValue + roadPartLength * ((roadRefVerts[i].z + roadRefLength / 2f) / roadRefLength)) / segLength;

                newMesh.vertsPos.Add(roadRefVerts[i]);
                newMesh.vertsPercLength.Add(vertsPercLength);
                newMesh.vertsPercWidth.Add(roadRefVerts[i].x / (roadRefWidth / 2f));                                  // -1 - 1
                newMesh.vertsPercHeight.Add((roadRefVerts[i].y + (roadRefHeight / 2f)) - roadRefHeight);
            }

            for (int i = 0; i < newMesh.tris.Count; i++)
            {
                newMesh.tris[i] += trisAddition;
            }

            return newMesh;
        }

        List<float> DistancedTs(int segIndex, List<float> prevTs) // prevTs is 0-1 in segment
        {
            List<float> newTs = new(prevTs);

            int rangeStart = indicesOfSegments[segIndex];
            int rangeEnd = segIndex < spline.NumSegments - 1 ? indicesOfSegments[segIndex + 1] : tsOfESP.Count;
            List<float> tsOfESPinSegment = tsOfESP.GetRange(rangeStart, rangeEnd - rangeStart);

            string output = "";
            prevTs.ForEach(x => output += x + ", ");
            Debug.Log(output);

            for (int i = 0; i < prevTs.Count; i++)
            {
                for (int t = 0; t < tsOfESPinSegment.Count; t++)
                {
                    if (prevTs[i] <= tsOfESPinSegment[t])
                    {
                        float prev_tsOfESPinSegment = t > 0 ? tsOfESPinSegment[t - 1] : 0;
                        newTs[i] = tsOfESPinSegment[t]; //(tsOfESPinSegment[t] - prev_tsOfESPinSegment) * 0.5f + prev_tsOfESPinSegment; //(prevTs[i] - prev_tsOfESPinSegment) / (tsOfESPinSegment[t] - prev_tsOfESPinSegment) + prev_tsOfESPinSegment;
                        break;
                    }
                }
            }

            output = "";
            newTs.ForEach(x => output += x + ", ");
            Debug.Log(output);

            return newTs;
        }

        #region Apply

        if (roadObject.GetComponent<MeshRenderer>() == null)
            roadObject.AddComponent<MeshRenderer>();
        if (roadObject.GetComponent<MeshFilter>() == null)
            roadObject.AddComponent<MeshFilter>();

        MeshFilter filter = roadObject.GetComponent<MeshFilter>();
        filter.sharedMesh = new Mesh() { vertices = mw.vertsPos.ToArray(), triangles = mw.tris.ToArray(), uv = mw.uvs.ToArray() };
        //GetComponent<MeshCollider>().sharedMesh = new Mesh() { vertices = mesh.vertsPos.ToArray(), triangles = mesh.tris.ToArray(), uv = mesh.uvs.ToArray() };
        filter.sharedMesh.RecalculateBounds();
        filter.sharedMesh.RecalculateNormals();

        roadObject.GetComponent<MeshRenderer>().materials = materials.ToArray();

        #endregion
    }

    class MeshWrapper
    {
        public List<Vector3> vertsPos;
        public List<float> vertsPercLength;
        public List<float> vertsPercWidth;
        public List<float> vertsPercHeight;
        public List<int> tris;
        public List<Vector2> uvs;

        public MeshWrapper()
        {
            vertsPos = new List<Vector3>();
            vertsPercHeight = new List<float>();
            vertsPercLength = new List<float>();
            vertsPercWidth = new List<float>();
            tris = new List<int>();
            uvs = new List<Vector2>();
        }

        public void Add(MeshWrapper newOne)
        {
            vertsPos.AddRange(newOne.vertsPos);
            vertsPercHeight.AddRange(newOne.vertsPercHeight);
            vertsPercWidth.AddRange(newOne.vertsPercWidth);
            vertsPercLength.AddRange(newOne.vertsPercLength);
            tris.AddRange(newOne.tris);
            uvs.AddRange(newOne.uvs);

            //Debug.Log("Added mw with " + newOne.vertsPos.Count + " vertices");
        }
    }

    */
}

public class MeshConfig
{
    public Mesh partMesh;
    public List<Material> materials;
    public float partLength;

    public MeshConfig(Mesh mesh, List<Material> materials, float partLength)
    {
        this.partMesh = mesh;
        this.materials = materials;
        this.partLength = partLength;
    }
}
