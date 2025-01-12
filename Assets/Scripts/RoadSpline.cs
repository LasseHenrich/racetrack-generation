using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

[Serializable]
public class RoadSpline : Spline
{
    public RoadSpline(Vector2 centre) : base(centre)
    {
    }

    public RoadSpline(List<Vector2> points, bool controlPointsIncluded = false, bool alreadyClosed = false) : base(points, controlPointsIncluded, alreadyClosed)
    {
    }

    #region Automatic Width Evaluation

    // Not used [16.02.2023]
    public float CalculateMinimumTPEradius()
    {
        var points = CalculateEvenlySpacedPoints(0.1f);

        float smallestRadius = float.MaxValue;

        for (int i = 0; i < points.Count; i++)
        {
            for (int j = 0; j < points.Count; j++)
            {
                if (Mathf.Min(Mathf.Abs(i - j + points.Count), Mathf.Abs(i - j - points.Count), Mathf.Abs(i - j)) < 5) continue;

                Vector2 disp = points[i] - points[j];
                Vector2 t_i = Tangent(i, points);

                Vector2 normal_proj = disp - Vector2.Dot(disp, t_i) * t_i;
                float numer = disp.sqrMagnitude;
                float denom = normal_proj.magnitude;
                float value = numer / denom;

                smallestRadius = Mathf.Min(smallestRadius, value);
            }
        }

        return smallestRadius;
    }

    // Used [16.02.2023]
    public float CalculateMinimumDistance()
    {
        float crossingPercRange = TopologyHandler.Crossing_EstimatedStartEndRange();
        crossingPercRange *= 2f; // Goodwill
        Debug.Log(crossingPercRange);

        var intersectingPercents = TopologyHandler.CalculateIntersectingPercents();
        var topologies = TopologyHandler.EmptyTopologyList();
        foreach (var ip in intersectingPercents)
        {
            var (a_start, a_end) = TopologyHandler.StartEnd01(ip.a - crossingPercRange * 0.5f, ip.a + crossingPercRange * 0.5f);
            var (b_start, b_end) = TopologyHandler.StartEnd01(ip.b - crossingPercRange * 0.5f, ip.b + crossingPercRange * 0.5f);

            TopologyHandler.AddTopologyToList(a_start, a_end, TopologyType.Crossing, ref topologies);
            TopologyHandler.AddTopologyToList(b_start, b_end, TopologyType.Crossing, ref topologies);
        }
        TopologyHandler.LogTopologies(topologies);


        float spacing = 1f;
        var points = CalculateEvenlySpacedPoints(spacing, out List<int> indicesOfSegments, out List<float> ts);

        float minPointOffset = 20 / spacing;
        List<int> indicesInCrossings = new();
        for (int i = 0; i < points.Count; i++)
            if (pointNearIntersection(i)) indicesInCrossings.Add(i);

        float shortestDist = float.MaxValue;

        //for (int i = 0; i < points.Count; i++)
        //{
        //    if (indicesInCrossings.Contains(i)) continue;
        //
        //    for (int j = i + 1; j < points.Count; j++) // j = i + 1: We don't want to calculate both (i, j) and (j, i)
        //    {
        //        if (Mathf.Min(Mathf.Abs(i - j + points.Count), Mathf.Abs(i - j - points.Count), Mathf.Abs(i - j)) < minPointOffset) continue;
        //        if (indicesInCrossings.Contains(j)) continue;
        //
        //
        //        float value = Vector2.Distance(points[i], points[j]);
        //
        //        //var a_tangent = Tangent(segmentIndexOfPoint(i), ts[i]);
        //        //var b_tangent = Tangent(segmentIndexOfPoint(j), ts[j]);
        //        //float tangentAngle = Vector2.Angle(a_tangent, b_tangent) * Mathf.Deg2Rad;
        //        //if (tangentAngle > Mathf.PI * 0.5) tangentAngle = Mathf.PI - tangentAngle;
        //        //float tanTangentAngle = Mathf.Tan(tangentAngle);
        //        //value *= 1 + tanTangentAngle;
        //
        //        //float value = radiusOfCircle(i, j);
        //
        //        shortestDist = Mathf.Min(shortestDist, value); 
        //    }
        //}

        //Debug.Log("shortestDist after point-to-point distance: " + shortestDist);

        var unitTangents = new List<Vector2>();
        for (int i = 0; i < points.Count; i++)
        {
            unitTangents.Add(UnitTangent(segmentIndexOfPoint(i), ts[i]));
        }
        var unitNormals = unitTangents.Select(t => t = new Vector2(t.y, -t.x)).ToList();

        minPointOffset = 5 / spacing;

        for (int i = 0; i < points.Count; i++)
        {
            if (indicesInCrossings.Contains(i)) continue;
            for (int j = i + 1; j < points.Count; j++)
            {
                if (indicesInCrossings.Contains(j)) continue;
                //if (Mathf.Min(Mathf.Abs(i - j + points.Count), Mathf.Abs(i - j - points.Count), Mathf.Abs(i - j)) < minPointOffset) continue;

                var (t, u) = MyMath.Intersection(points[i], unitNormals[i], points[j], unitNormals[j]);

                t = Mathf.Abs(t); u = Mathf.Abs(u);
                float value = (t + u) * 0.5f * 2f; // Mean of t and u, but *2 as we need the whole distance

                //if (t < shortestDist)
                //    Debug.Log("New shortest Dist: t=" + t + " from " + points[i] + " to " + points[j] + " where i=" + i + " and j=" + j + ". Tangents are " + tangents[i] + " and " + tangents[j] + ". Normals are " + normals[i] + " and " + normals[j]);
                //if (u < shortestDist)
                //    Debug.Log("New shortest Dist: u=" + u + " from " + points[i] + " to " + points[j] + " where i=" + i + " and j=" + j + ". Tangents are " + tangents[i] + " and " + tangents[j] + ". Normals are " + normals[i] + " and " + normals[j]);

                //if (value < shortestDist)
                //    Debug.Log("New shortest Dist: value=" + value + " from " + points[i] + " to " + points[j] + " where i=" + i + " and j=" + j + ". Tangents are " + tangents[i] + " and " + tangents[j] + ". Normals are " + normals[i] + " and " + normals[j]);

                //float tangentAngle = Vector2.Angle(unitTangents[i], unitTangents[j]);
                //float tanTangentAngle = Mathf.Tan(tangentAngle);
                //value *= 1 + tangentAngle;

                shortestDist = Mathf.Min(shortestDist, value);
            }
        }

        //Debug.Log("shortestDist after normal-normal intersection: " + shortestDist);

        int segmentIndexOfPoint(int p_index)
        {
            for (int segIndex = 1; segIndex < indicesOfSegments.Count; segIndex++)
            {
                if (p_index < indicesOfSegments[segIndex])
                    return segIndex - 1;
            }

            return indicesOfSegments.Count - 1;
        }

        bool pointNearIntersection(int p_index)
        {
            float dstAlongSpline = GetDstAlongSpline(segmentIndexOfPoint(p_index), ts[p_index], false);
            return TopologyHandler.GetTopologyAndTInTopology(dstAlongSpline, topologies).type == TopologyType.Crossing;
        }

        /*

        float radiusOfCircle(int a, int b)
        {
            var a_pt = points[a];
            var b_pt = points[b];

            Vector2 disp = a_pt - b_pt;
            Vector2 a_tangent = Tangent(segmentIndexOfPoint(a), ts[a]);

            Vector2 normal_proj = disp - Vector2.Dot(disp, a_tangent) * a_tangent;
            float numer = Mathf.Pow(disp.magnitude, 2);
            float denom = Mathf.Abs(normal_proj.magnitude);
            return numer / denom;
        }
        */

        return shortestDist;
    }

    #endregion

    #region Creating the Mesh

    // Invserse Function of GetDstAlongSpline
    public (Vector2 point, int segIndex, float t) CalculatePointAtDist(float dist, float resolution = 1)
    {
        if (dist >= TotalLength)
        {
            if (isClosed) dist %= TotalLength;
            else return GetLastPoint();
        }
        else if (dist < 0)
        {
            if (isClosed) dist = (dist % TotalLength + TotalLength);
            else return GetFirstPoint();
        }

        int segIndex;
        float prevSegsLength = 0;
        for (segIndex = 0; segIndex < NumSegments; segIndex++)
        {
            float nextSegLength = GetSegmentLength(segIndex);
            if (prevSegsLength + nextSegLength > dist)
            {
                break;
            }
            prevSegsLength += nextSegLength;
        }

        int divisions = (int)(10 * GetSegmentLength(segIndex) * resolution);
        Vector2[] p = GetPointsInSegment(segIndex);
        Vector2 previousPoint = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], 0);
        float totalMeasuredLength = prevSegsLength;
        float t = 0;
        Vector2 point = p[0];

        if (totalMeasuredLength < dist) // Otherwise, t = 0
        {
            for (int i = 0; i < divisions; i++)
            {
                t = (i + 1) / (float)divisions;
                point = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], t);
                totalMeasuredLength += Vector2.Distance(point, previousPoint);
                previousPoint = point;
                if (totalMeasuredLength > dist)
                    break;
            }
        }

        return (point, segIndex, t);


        //return CalculatePointAfterRef(spacingToRef: dist, 0, 0, resolution);
    }

    public (Vector2 point, int segIndex, float t) GetLastPoint()
    {
        if (IsClosed)
            return GetFirstPoint(); // For more accurate looping. I don't know why, but this is slightly different
        return (controls[^1], NumSegments - 1, 1);
    }

    public (Vector2 point, int segIndex, float t) GetFirstPoint()
    {
        return (FirstPoint, 0, 0);
    }

    /// <summary>
    /// Calculates and returns a point, that is spacingToPrev on the curve after refPoint
    /// </summary>
    /// <param name="spacingToRef"></param>
    /// <param name="refSegIndex">Index of segment of previous point</param>
    /// <param name="refT">t of previous point</param>
    /// <param name="resolution"></param>
    /// <returns>point, segmentIndex, t</returns>
    public (Vector2 point, int segIndex, float t) CalculatePointAfterRef(float spacingToRef, int refSegIndex, float refT, float resolution = 1)
    {
        // Note: This algorithm should be adapted to not calculate through every possible sement until the point
        // Instead, we can just skip all segments that are clearly before the point and then only calculate within the next segment

        if (spacingToRef > TotalLength)
            spacingToRef -= TotalLength;

        Vector2[] p_ref = GetPointsInSegment(refSegIndex);
        Vector2 refPoint = Bezier.EvaluateCubic(p_ref[0], p_ref[1], p_ref[2], p_ref[3], refT);

        float t = refT;
        bool looped = false;
        float spacing = 0;
        Vector2 prev_pointOnCurve = refPoint;

        for (int segmentIndex = refSegIndex; segmentIndex < NumSegments; segmentIndex++)
        {
            Vector2[] p = GetPointsInSegment(segmentIndex);
            float segmentLength = GetSegmentLength(segmentIndex);
            int divisions = Mathf.CeilToInt(segmentLength * resolution * 10);

            while (t <= 1f)
            {

                Vector2 pointOnCurve = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], t);
                float dstTo_prev_pointOnCurve = Vector2.Distance(prev_pointOnCurve, pointOnCurve);
                spacing += dstTo_prev_pointOnCurve;

                if (spacing >= spacingToRef)
                {
                    Vector2 point = pointOnCurve;
                    float overshootDst = spacing - spacingToRef;
                    point -= (pointOnCurve - refPoint).normalized * overshootDst;
                    if (spacing > 0)
                    {
                        float overshootT = (overshootDst / dstTo_prev_pointOnCurve) / divisions; // Seems to be more precise than overshootDst / segmentLength. With the latter one we have very small gaps between segments
                        t -= overshootT;
                    }

                    return (point, segmentIndex, t);
                }

                prev_pointOnCurve = pointOnCurve;
                t += 1f / divisions;
            }

            t = 0;

            if (segmentIndex == NumSegments - 1)
            {
                if (looped)
                {
                    Debug.LogError("Something went wrong: Looped two times during point calculation. Maybe spacingToPrev is too large?");
                    return (Vector2.zero, -1, -1);
                }
                looped = true;
                segmentIndex = -1; // (The loop will set it to 0)
            }
        }

        Debug.LogError("Something went wrong: Shouldn't be here.");
        return (Vector2.zero, -1, -1);
    }

    #endregion

    #region Intersections

    /// <summary>
    /// Returns a new RoadSpline that's a copy of the current but with an intersection (more).
    /// </summary>
    /// <param name="preferredDistance">Defines how far long the straight section around a crossing should be. Usually, 0 is a good value.</param>
    /// <returns></returns>
    public void AddIntersections(float preferredDistance)
    {
        float spacing = 0.2f;
        List<Vector2> points = CalculateEvenlySpacedPoints(spacing, out List<int> indicesOfSegments);

        List<CrossingCandidate> candidates = AssembleCrossingCandidates(points, spacing);

        if (candidates.Count == 0)
            Debug.Log("Didn't find any candidates.");
        else Debug.Log("Number of candidates: " + candidates.Count);

        //string output = "";
        //candidates.ForEach(c => output += (points[c.i.index] + " to " + points[c.j.index]) + "\n");
        //Debug.Log(output);

        List<CrossingCandidatePair> candidatePairs = AssembleCandidatePairs(preferredDistance, points, indicesOfSegments, candidates);

        var newRoadSpline = this;

        if (candidatePairs.Count > 0)
        {
            var ccp = candidatePairs.OrderBy(x => x.badness).First();
            //Debug.Log("Min badness: " + ccp.badness);
            controls = new(this, ccp.spline.controls.GetList());
            CalculateAllPointsInSegments();
            CalculateAllSegmentLengths(); // THOSE ARE NECESSARY, the one in the constructor doesn't seem to work [2023-03-26]
        }
        else
        {
            Debug.LogWarning("No candidate pairs found.");
        }
    }

    private List<CrossingCandidatePair> AssembleCandidatePairs(float preferredDistance, List<Vector2> points, List<int> indicesOfSegments, List<CrossingCandidate> candidates)
    {
        List<CrossingCandidatePair> candidatePairs = new();

        for (int c1 = 0; c1 < candidates.Count; c1++)
        {
            for (int c2 = c1 + 1; c2 < candidates.Count; c2++)
            {
                var cc1 = candidates[c1];
                var cc2 = candidates[c2];

                //Debug.Log("Testing New pair: " + cc1.i + ", " + cc1.j + ", " + cc2.i + ", " + cc2.j
                //    + " with points " + points[cc1.i.index] + ", " + points[cc1.j.index] + ", " + points[cc2.i.index] + ", " + points[cc2.j.index]);

                List<CrossingOpeningPoint> sortCOPs = new() { cc1.i, cc1.j, cc2.i, cc2.j }; // Sorted CrossingOpeningPoints
                sortCOPs.Sort();

                if (!((sortCOPs[0] == cc1.i || sortCOPs[0] == cc1.j) && (sortCOPs[3] == cc1.i || sortCOPs[3] == cc1.j)) &&    // not 1221
                    !((sortCOPs[0] == cc2.i || sortCOPs[0] == cc2.j) && (sortCOPs[3] == cc2.i || sortCOPs[3] == cc2.j)) &&    // not 2112
                    !((sortCOPs[0] == cc1.i || sortCOPs[0] == cc1.j) && (sortCOPs[1] == cc1.i || sortCOPs[1] == cc1.j)) &&    // not 1122
                    !((sortCOPs[0] == cc2.i || sortCOPs[0] == cc2.j) && (sortCOPs[1] == cc2.i || sortCOPs[1] == cc2.j)))      // not 2211
                {
                    //Debug.Log("Pair-elements are ordered the right way (1212 or 2121)");

                    if (sortCOPs[0].forward != sortCOPs[1].forward && sortCOPs[0].forward != sortCOPs[2].forward && sortCOPs[0].forward == sortCOPs[3].forward)  // fbbf or bffb
                                                                                                                                                                 //|| (sortCOPs[0].forward == sortCOPs[1].forward && sortCOPs[0].forward != sortCOPs[2].forward && sortCOPs[0].forward != sortCOPs[3].forward))   // ffbb or bbff
                    {
                        //Debug.Log("Pair-elementts are matchable (fbbf or bffb");//, ffbb or bbff)");

                        var (c1_t, c2_t) = MyMath.Intersection(points, cc1.i.index, cc1.j.index, cc2.i.index, cc2.j.index); // Use different method signature, as tangent my not be pointing towards other point but rather away from it

                        // Do segments overlap?
                        if (c1_t >= 0 && c1_t <= 1 && c2_t >= 0 && c2_t <= 1)
                        {
                            //Debug.Log("Pair-elements overlap");

                            RoadSpline splineCandidate = GenerateIntersectingSplineCandidate(points, indicesOfSegments, sortCOPs);

                            // How they overlap -> how well would a crossing be
                            (bool acceptable, float badness) = QualifyIntersectingSplineCandidate(preferredDistance, points, cc1, cc2, splineCandidate);

                            //Debug.Log("Badness: " + badness);

                            candidatePairs.Add(new(badness, splineCandidate));
                        }
                    }
                }
            }
        }

        return candidatePairs;
    }

    /// <returns><acceptable, badness</returns>
    private Tuple<bool, float> QualifyIntersectingSplineCandidate(float preferredDistance, List<Vector2> points, CrossingCandidate cc1, CrossingCandidate cc2, RoadSpline splineCandidate)
    {
        // ToDo: Return false as first element if outside of certain bounds

        float badness = 0;
        float sqrPreferredDistance = preferredDistance * preferredDistance;
        badness += Mathf.Abs(Vector2.SqrMagnitude(points[cc1.i.index] - points[cc1.j.index]) - sqrPreferredDistance); // Impact of distance of cc1.i and cc1.j -> The closer the better
        badness += Mathf.Abs(Vector2.SqrMagnitude(points[cc2.i.index] - points[cc2.j.index]) - sqrPreferredDistance); // Impact of distance of cc2.i and cc2.j -> The closer the better
        badness += Mathf.Abs(TotalLength - splineCandidate.TotalLength); // Length of curve should not become much smaller -> We don't wanna lose much information / turns
        var a_tangent = points[cc1.i.index] - points[cc1.j.index];
        var b_tangent = points[cc2.i.index] - points[cc2.j.index];
        float angle = Vector2.Angle(a_tangent, b_tangent) * Mathf.Deg2Rad;  // [0, Mathf.PI)
        if (angle > Mathf.PI * 0.5f)
            angle = Mathf.PI - angle;                                       // [0, Mathf.PI * 0.5f)
        float invAngle = Mathf.PI * 0.5f - angle;                           // Inversed
        float tanInvAngle = Mathf.Tan(invAngle);
        badness += tanInvAngle * 50;

        return Tuple.Create(true, badness);
    }

    private RoadSpline GenerateIntersectingSplineCandidate(List<Vector2> points, List<int> indicesOfSegments, List<CrossingOpeningPoint> sortCOPs)
    {
        var splineCandidate = DeepCopy();
        List<int> c_indicesOfSegments = new(indicesOfSegments);

        // 1. Split Segments at cops, so that cops are anchors
        List<int> copAnchorIndices = new();
        int copIndex = 0;
        CrossingOpeningPoint cop = sortCOPs[copIndex];

        for (int segIndex = 0; segIndex < splineCandidate.NumSegments; segIndex++)
        {
            // Is point in previous segment?
            if (segIndex == splineCandidate.NumSegments - 1 || cop.index < c_indicesOfSegments[segIndex + 1]) // Subtract copAnchorIndices.Count, as we added copAnchorIndices.Count new segments with SplitSegment
            {
                int pointsInSeg = (segIndex < splineCandidate.NumSegments - 1 ? c_indicesOfSegments[segIndex + 1] : points.Count) - c_indicesOfSegments[segIndex];
                float tInSeg = (cop.index - c_indicesOfSegments[segIndex]) / (float)pointsInSeg;
                splineCandidate.SplitSegment(points[cop.index], segIndex, tInSeg: tInSeg); // We have copAnchorIndices.Count segments more through previous splits
                                                                                           //Debug.Log("new cop " + points[cop.index] + " in seg " + segIndex);
                                                                                           //Debug.Log("tInSeg: " + tInSeg);
                copAnchorIndices.Add(segIndex * 3 + 3);

                c_indicesOfSegments.Insert(segIndex + 1, cop.index);
                copIndex++;
                if (copIndex == sortCOPs.Count)
                    break;
                cop = sortCOPs[copIndex];
            }
        }

        // 2. Create new control-List
        // Note: Origin is points[0], End is points[^1]. If you " + 1 + 1" it means: +1 for Index->Count and +1 for next control point
        List<Vector2> newControls = new();
        if (sortCOPs[0].forward != sortCOPs[1].forward && sortCOPs[0].forward != sortCOPs[2].forward && sortCOPs[0].forward == sortCOPs[3].forward)  // fbbf or bffb
        {
            if (!sortCOPs[0].forward && sortCOPs[1].forward && sortCOPs[2].forward && !sortCOPs[3].forward) // bffb
            {
                /*
                 * 1 to 0, 0 to 2, 2 to 3, 3 to 1
                 */
                //Debug.Log("bffb");

                Tuple<int, int> cop1To0 = new(
                    copAnchorIndices[0],
                    (copAnchorIndices[1] - copAnchorIndices[0]) + 1 // Must start with anchor
                );
                AddControlRange(cop1To0, true);

                Tuple<int, int> copBefore0 = new(
                    LoopIndex(copAnchorIndices[0] - 1),
                    1
                );
                AddControlRange(copBefore0, false);

                Tuple<int, int> cop2To3 = new(
                    copAnchorIndices[2] - 1,
                    (copAnchorIndices[3] - (copAnchorIndices[2] - 1)) + 1 + 1
                );
                AddControlRange(cop2To3, false);

                Tuple<int, int> copAfter1 = new(
                    copAnchorIndices[1] + 1,
                    1
                );
                AddControlRange(copAfter1, false);
            }
            else if (sortCOPs[0].forward && !sortCOPs[1].forward && !sortCOPs[2].forward && sortCOPs[3].forward) // fbbf
            {
                /*
                 * 0 to Origin, Origin to End, End to 3, 3 to 1, 1 to 2, 2 to 0
                 */
                //Debug.Log("fbbf");

                Tuple<int, int> cop0ToOrigin = new(
                    0,
                    (copAnchorIndices[0] - 0) + 1
                );
                AddControlRange(cop0ToOrigin, true);

                Tuple<int, int> copEndTo3 = new(
                   copAnchorIndices[3] - 1,
                   ((splineCandidate.controls.Count - 1) - (copAnchorIndices[3] - 1)) + 1
                );
                AddControlRange(copEndTo3, true);

                Tuple<int, int> cop1To2 = new(
                    copAnchorIndices[1] - 1,
                    (copAnchorIndices[2] - (copAnchorIndices[1] - 1)) + 1 + 1
                );
                AddControlRange(cop1To2, false);

                Tuple<int, int> copAfter0 = new(
                    copAnchorIndices[0] + 1,
                    1
                );
                AddControlRange(copAfter0, false);
            }
            else
            {
                Debug.LogWarning("Case not yet covered: " + sortCOPs[0].forward + " " + sortCOPs[1].forward + " " + sortCOPs[2].forward + " " + sortCOPs[3].forward);
            }
        }
        else // ffbb or bbff
        {
            Debug.Log("Case not yet covered. (Couldn't find a suitable sketch.)");
            // Edit: after lots of track generations this hasn't happened yet, so likely impossible.
        }

        void AddControlRange(Tuple<int, int> startCount, bool reversed)
        {
            List<Vector2> toAdd = splineCandidate.controls.GetRange(startCount.Item1, startCount.Item2);
            if (reversed) toAdd.Reverse();
            newControls.AddRange(toAdd);
        }

        splineCandidate.controls.RemoveRange(0, splineCandidate.controls.Count);
        splineCandidate.controls.AddRange(newControls);
        return splineCandidate;
    }

    private List<CrossingCandidate> AssembleCrossingCandidates(List<Vector2> points, float spacing)
    {
        List<CrossingCandidate> candidates = new(); // <i, j, i_forward, j_forward>

        List<Vector2> tangents = GetTangents(points);
        int minIndicesDist = (int)(10 / spacing);                // Minimum Distance betweeen Indices of points of one pair
        float maxTangentCrossThreshold = 0.05f; // Every Cross(tan[i],-tan[j]) below this number will be a valid candidate
        float maxDispCrossTheshold = 0.08f;     // Every Cross(unitDisp, (tan[i]+tan[j]).norm) below this number will be a valid candidate
        int skipIncidesAfterFound = 1;//(int)(1f / spacing);         // Indices to be skipped after a valid candidate was found. Note that I've also experimented with making this variable based on the change in tangent from i to i + 1. But, as spacing to small, this doesn't really have any effect.


        for (int i = 0; i < points.Count; i++)
        {
            for (int j = i + 1; j < points.Count; j++)
            {
                // Are points far enough apart along the curve?
                if (j - i >= minIndicesDist)
                {
                    //Debug.Log("Points " + i + " and " + j + " are far enough apart: " + points[i] + " & " + points[j]);

                    float tangentCross = Mathf.Abs(CurveGenUtils.Cross(tangents[i], -tangents[j]));

                    // Are lines roughly parallel?
                    if (tangentCross < maxTangentCrossThreshold)
                    {
                        //Debug.Log("Tangents rougly parallel: " + tangents[i] + " & " + tangents[j]);

                        Vector2 unitDisp_i2j = (points[j] - points[i]).normalized;
                        Vector2 meanTangent = (tangents[i] + tangents[j]); // At this point, we don't know whether it points from i to j or j to i
                        if (meanTangent.sqrMagnitude < 0.5f) meanTangent = tangents[i] - tangents[j]; // Should take the inverse of one tangent, because they point in exactly the opposite direction
                        meanTangent.Normalize();

                        if (meanTangent.sqrMagnitude < 0.9f)
                        {
                            Debug.LogWarning("Warning: Got meanTangent.magnitude < 1: " + meanTangent.sqrMagnitude + ". Skipping...");
                            continue;
                        }

                        float dispCross = Mathf.Abs(CurveGenUtils.Cross(unitDisp_i2j, meanTangent));

                        // Do tangents roughly point at each other? Taking the position of the points into account
                        if (dispCross < maxDispCrossTheshold)
                        {
                            //Debug.Log("Points roughly point at each other. Valid candidate found: " + points[i] + " & " + points[j] + ". Tangents: " + tangents[i] + " & " + tangents[j] + " meanTangent was " + meanTangent);

                            if (i == 154 && j == 278)
                            {
                                Debug.Log(meanTangent);
                            }

                            bool i_forward = Vector2.Dot(unitDisp_i2j, (points[LoopIndex(i + 1, points)] - points[i]).normalized) > 0; // tangent points in direction of next point
                            bool j_forward = Vector2.Dot(unitDisp_i2j, (points[LoopIndex(j + 1, points)] - points[j]).normalized) > 0;

                            // Check if both tangents point in right direction on curve, i.e. forward or backwards
                            candidates.Add(new() { i = new() { index = i, forward = i_forward }, j = new() { index = j, forward = j_forward } });

                            // If we found a candidate, chances are high we'll find another one with the next j, so we skip a few
                            // Also we don't want any too near candidates, so we have to skip because of this as well
                            i += skipIncidesAfterFound;
                            if (i > points.Count - 1)
                                break; // break is sufficient, as i-loop will then break again
                        }
                    }
                }
            }
        }

        return candidates;
    }

    /// <summary>
    /// Better version. See https://pomax.github.io/bezierinfo/#splitting.
    /// </summary>
    public void SplitSegment(Vector2 anchorPos, int segmentIndex, float tInSeg)
    {
        int prevAnchorPoint = segmentIndex * 3;

        Vector2 midBeforeCtrlPts = // Mid of two oiginal control points (between green and blue)
            Between(prevAnchorPoint + 1, prevAnchorPoint + 2);

        controls[LoopIndex(prevAnchorPoint + 1)] = // set ctrl point after prevAnchorPoint (green to between red and green)
            Between(prevAnchorPoint, prevAnchorPoint + 1);
        controls[LoopIndex(prevAnchorPoint + 2)] = // set ctrl point before enxt anchor point (blue to between yellow and blue)
            Between(prevAnchorPoint + 2, prevAnchorPoint + 3);

        controls.InsertRange(segmentIndex * 3 + 2, new Vector2[] {
            BetweenV(controls[LoopIndex(prevAnchorPoint + 1)], midBeforeCtrlPts), // Set first ctrl of new segment
            anchorPos,
            BetweenV(midBeforeCtrlPts, controls[LoopIndex(prevAnchorPoint + 2)])
        });

        Vector2 Between(int a, int b)
        {
            return BetweenV(controls[LoopIndex(a)], controls[LoopIndex(b)]);
        }

        Vector2 BetweenV(Vector2 A, Vector2 B)
        {
            return A + (B - A) * tInSeg;
        }

        //segmentTypes.Insert(segmentIndex, segmentTypes[segmentIndex]); // Should also work without [2023/03/09]
    }

    /// <summary>
    /// Better version. See https://pomax.github.io/bezierinfo/#splitting. Returns two segments à 4 points.
    /// </summary>
    public (List<Vector2>, List<Vector2>) SplitSegment(List<Vector2> controls, float tInSeg)
    {
        List<Vector2> newControls = new()
        {
            controls[0],
            Between(0, 1),
            BetweenV(Between(0, 1), Between(1, 2)),
            BetweenV(BetweenV(Between(0, 1), Between(1, 2)), BetweenV(Between(1, 2), Between(2, 3))),
            BetweenV(Between(1, 2), Between(2, 3)),
            Between(2, 3),
            controls[3]
        };

        Vector2 Between(int a, int b)
        {
            return BetweenV(controls[a], controls[b]);
        }

        Vector2 BetweenV(Vector2 A, Vector2 B)
        {
            return Vector2.Lerp(A, B, tInSeg);
        }

        return (newControls.GetRange(0, 4), newControls.GetRange(3, 4));
    }

    #endregion

    #region Features

    /*
    public void AddFeatures(float minBridgeStartEndLength)
    {
        var intersections = CalculateIntersections();
        intersections.ForEach(x => Debug.Log(x));

        
        segmentTypes = Enumerable.Repeat(SegmentType.Normal, NumSegments).ToList();

        foreach (var intersection in intersections)
        {
            var a = intersection.Item1;
            var b = intersection.Item2;

            if (IntersectionsNormal(a, b))
            {
                var ran = UnityEngine.Random.value;
                if (ran < 0.5f)
                    AddShortBridge(a);
                else
                    AddShortBridge(b);
            }
            else
            {
                Debug.LogWarning("This will be interesting... Do not override curve Lasse!");
            }
        }

        bool IntersectionsNormal(int a, int b)
        {
            return segmentTypes[LoopIndex(a - 1, segmentTypes)] == SegmentType.Normal && segmentTypes[a] == SegmentType.Normal && segmentTypes[LoopIndex(a + 1, segmentTypes)] == SegmentType.Normal &&
                   segmentTypes[LoopIndex(b - 1, segmentTypes)] == SegmentType.Normal && segmentTypes[b] == SegmentType.Normal && segmentTypes[LoopIndex(b + 1, segmentTypes)] == SegmentType.Normal;
        }

        void AddShortBridge(int centerSeg)
        {
            int bridgeStartSeg = BestBridgeStartEnd(-1, -5);
            int bridgeEndSeg = BestBridgeStartEnd(+1, 5);

            segmentTypes[bridgeStartSeg] = SegmentType.BridgeStart;
            segmentTypes[bridgeEndSeg] = SegmentType.BridgeEnd;

            SegmentsInBetween(bridgeStartSeg, bridgeEndSeg).ForEach(
                x => segmentTypes[x] = SegmentType.Bridge
            );



            int BestBridgeStartEnd(int increment, int maxIncrement)
            {
                List<Tuple<int, float>> searchedSegments = new();
                int currIncrement = increment;
                while (currIncrement != maxIncrement)
                {
                    int segIndex = LoopIndex(centerSeg + currIncrement, NumSegments);
                    float length = SegmentLength(segIndex);

                    if (length >= minBridgeStartEndLength)
                        return segIndex;

                    searchedSegments.Add(new Tuple<int, float>(segIndex, length));
                    currIncrement += increment;
                }

                return searchedSegments.OrderBy(x => x.Item2).Last().Item1;
            }
        }
        
    }
    */

    public List<Tuple<Tuple<int, float>, Tuple<int, float>>> CalculateIntersections()
    {
        const float threshold = 0.001f;
        List<Tuple<Tuple<int, float>, Tuple<int, float>>> intersections = new();

        for (int s1_index = 0; s1_index < NumSegments; s1_index++)
        {
            for (int s2_index = s1_index + 1; s2_index < NumSegments; s2_index++)
            {
                var s1_initialPts = GetPointsInSegment(s1_index).ToList();
                var s2_initialPts = GetPointsInSegment(s2_index).ToList();

                var s1_initialStartEnd = new Tuple<float, float>(0f, 1f);
                var s2_initialStartEnd = new Tuple<float, float>(0f, 1f);

                if (s2_index == s1_index + 1)
                {
                    (s1_initialPts, _) = SplitSegment(s1_initialPts, 0.9f);
                    s1_initialStartEnd = new(0f, 0.9f);
                }
                else if (s1_index == 0 && s2_index == NumSegments - 1)
                {
                    (_, s1_initialPts) = SplitSegment(s1_initialPts, 0.1f);
                    s2_initialStartEnd = new(0.1f, 1f);
                }

                var (intersecting, s1, s2)
                    = BezInt(new(s1_initialPts, s1_initialStartEnd), new(s2_initialPts, s2_initialStartEnd));

                if (intersecting)
                {
                    float s1_t = s1.startEndOfOrigSeg.Item1 + (s1.startEndOfOrigSeg.Item2 - s1.startEndOfOrigSeg.Item1) * 0.5f;
                    float s2_t = s2.startEndOfOrigSeg.Item1 + (s2.startEndOfOrigSeg.Item2 - s2.startEndOfOrigSeg.Item1) * 0.5f;
                    intersections.Add(new(new(s1_index, s1_t), new(s2_index, s2_t)));

                    //Vector2 pos2D = CalculatePointAtT(s1_index, s1_t);
                    //Vector3 pos3D = new(pos2D.x, 0, pos2D.y);
                    //GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                    //cube.transform.position = new(s1.ctrlPts[0].x, 0, s1.ctrlPts[0].y);
                    //
                    //pos2D = CalculatePointAtT(s1_index, s1.startEndOfOrigSeg.Item1);
                    //pos3D = new(pos2D.x, 0, pos2D.y);
                    //cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                    //cube.transform.position = pos3D;

                    //List<Vector2> pts = s1_initialPts;
                    //List<Vector2> pts2;
                    //for (int i = 0; i < 1; i++)
                    //{
                    //    (pts, pts2) = SplitSegment(pts, 0.5f);
                    //    pts.ForEach(pt => Place(pt));
                    //    pts2.ForEach(pt => Place(pt));
                    //}

                    //var (a_ctrl, _) = SplitSegment(s1_initialPts, 0.5f);
                    //Place(a_ctrl[3]);
                    //Place(CalculatePointAtT(s1_index, 0.5f));

                    //void Place(Vector2 pt)
                    //{
                    //    GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                    //    cube.transform.position = new(pt.x, 0, pt.y);
                    //}
                }
            }
        }

        // https://stackoverflow.com/questions/4039229/checking-if-two-cubic-b%C3%A9zier-curves-intersect
        (bool, SubSegment, SubSegment) BezInt(SubSegment s1, SubSegment s2)
        {
            var (bbox1_min, bbox1_max) = BoundingBox(s1.ctrlPts);
            var (bbox2_min, bbox2_max) = BoundingBox(s2.ctrlPts);

            bool bboxesOverlapping = bbox1_max.x >= bbox2_min.x && bbox2_max.x >= bbox1_min.x &&
                                     bbox1_max.y >= bbox2_min.y && bbox2_max.y >= bbox1_min.y;
            if (!bboxesOverlapping)
                return (false, null, null);

            var area_bbox1 = Vector2.SqrMagnitude(bbox1_max - bbox1_min);
            var area_bbox2 = Vector2.SqrMagnitude(bbox2_max - bbox2_min);
            if (area_bbox1 + area_bbox2 < threshold)
                return (true, s1, s2);

            var (s1a_ctrl, s1b_ctrl) = SplitSegment(s1.ctrlPts, 0.5f);
            var (s2a_ctrl, s2b_ctrl) = SplitSegment(s2.ctrlPts, 0.5f);

            var (s1a, s1b) = s1.SplitAtMid(this);
            var (s2a, s2b) = s2.SplitAtMid(this);

            var pairs = new List<Tuple<SubSegment, SubSegment>>() { new(s1a, s2a), new(s1a, s2b), new(s1b, s2a), new(s1b, s2b) };

            foreach (var pair in pairs)
            {
                var (intersecting, new_s1, new_s2) = BezInt(pair.Item1, pair.Item2);
                if (intersecting)
                {
                    return (true, new_s1, new_s2);
                }
            }
            return (false, null, null);
        }

        return intersections;
    }

    class SubSegment
    {
        public List<Vector2> ctrlPts;
        public Tuple<float, float> startEndOfOrigSeg;

        public SubSegment(List<Vector2> ctrlPts, Tuple<float, float> startEndOfOrigSeg)
        {
            this.ctrlPts = ctrlPts;
            this.startEndOfOrigSeg = startEndOfOrigSeg;
        }

        public (SubSegment, SubSegment) SplitAtMid(RoadSpline spline)
        {
            var (a_ctrl, b_ctrl) = spline.SplitSegment(ctrlPts, 0.5f);

            var start = startEndOfOrigSeg.Item1;
            var mid = startEndOfOrigSeg.Item1 + (startEndOfOrigSeg.Item2 - startEndOfOrigSeg.Item1) * 0.5f;
            var end = startEndOfOrigSeg.Item2;
            var a_startEnd = new Tuple<float, float>(start, mid);
            var b_startEnd = new Tuple<float, float>(mid, end);

            //Debug.Log(startEndOfOrigSeg + " -> " + a_startEnd + " and " + b_startEnd);

            return (new(a_ctrl, a_startEnd), new(b_ctrl, b_startEnd));
        }
    }

    /*
    public SegmentType GetSegmentType(int segIndex)
    {
        if (segmentTypes == null)
        {
            //Debug.LogWarning("segmentTypes not initialized... returning SegmentType.Normal.");
            return SegmentType.Normal;
        }
        return segmentTypes[segIndex];
    }
    */

    #endregion

    #region Helper Functions

    public Vector2 UnitTangent(int segmentIndex, float t)
    {
        /*
        float tDiff = 0.01f;

        Vector2[] p = GetPointsInSegment(segmentIndex);
        Vector2 point = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], t);

        int prevSegmentIndex = segmentIndex;
        float tPrev = t - tDiff;
        if (tPrev < 0)
        {
            tPrev += 1;
            prevSegmentIndex = (prevSegmentIndex - 1 + NumSegments) % NumSegments;
        }
        p = GetPointsInSegment(prevSegmentIndex);
        Vector2 prev = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], tPrev);

        int nextSegmentIndex = segmentIndex;
        float tNext = t + tDiff;
        if (tNext > 1)
        {
            tNext -= 1;
            nextSegmentIndex = (nextSegmentIndex + 1) % NumSegments;
        }
        p = GetPointsInSegment(nextSegmentIndex);
        Vector2 next = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], tNext);

        return ((point - prev) + (next - point)).normalized;
        */

        Vector2[] p = GetPointsInSegment(segmentIndex);
        Vector3 tangent3D = Bezier.Tangent(p[0], p[1], p[2], p[3], t);
        return new Vector2(tangent3D.x, tangent3D.y).normalized;
    }

    Vector2 Tangent(int i, List<Vector2> points)
    {
        Vector2 prev = points[(i - 1 + points.Count) % points.Count];
        Vector2 next = points[(i + 1) % points.Count];
        return ((points[i] - prev) + (next - points[i])).normalized;
    }

    List<Vector2> GetTangents(List<Vector2> points)
    {
        return points.Select((p, i) => Tangent(i, points)).ToList();
    }

    public Vector2 FirstPoint
    {
        get
        {
            return controls[0];
        }
    }

    private List<int> SegmentsInBetween(int from, int to)
    {
        List<int> segments = new();

        int currIndex = from;
        currIndex = LoopIndex(currIndex + 1, NumSegments);
        while (currIndex != to)
        {
            segments.Add(currIndex);
            currIndex = LoopIndex(currIndex + 1, NumSegments);
        }

        return segments;
    }

    public Vector2 Normal(int segIndex, float t)
    {
        var tangent = UnitTangent(segIndex, t);
        var normal = new Vector2(-tangent.y, tangent.x);
        return normal;
    }

    public RoadSpline DeepCopy()
    {
        var copy = new RoadSpline(new(controls.GetList()), true, IsClosed);
        return copy;
    }

    #endregion
}

class CrossingOpeningPoint : IComparable<CrossingOpeningPoint>
{
    public int index;
    public bool forward;

    public int CompareTo(CrossingOpeningPoint other)
    {
        return index.CompareTo(other.index);
    }

    public override string ToString()
    {
        return index.ToString();
    }
}

class CrossingCandidate
{
    public CrossingOpeningPoint i, j;

    public float SqrSegLength(List<Vector2> points)
    {
        return Vector2.SqrMagnitude(points[i.index] - points[j.index]);
    }
}

class CrossingCandidatePair
{
    public float badness;
    public RoadSpline spline;

    public CrossingCandidatePair(float badness, RoadSpline spline)
    {
        this.badness = badness;
        this.spline = spline;
    }
}