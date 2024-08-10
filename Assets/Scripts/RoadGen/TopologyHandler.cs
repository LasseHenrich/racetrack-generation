using MathNet.Numerics.Differentiation;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Net.WebSockets;
using UnityEngine;

public static class TopologyHandler
{
    public static List<TopologyRange> topologies;
    public static List<List<float>> fourPercCrossings; // A list of crossings. A crossing is represented by four percs, where their order is equal to the visual order of the real quandrants, CLOCKWISE!

    private static List<(float start, float end)> intersectingPercents;

    private static float splineLength;

    static float roadPartLength;
    static float roadPartWidth;
    static float bridgeAscentLength;
    static float bridgeDescentLength;
    static float crossingExtraSize;
    static float rampLength;
    static RoadSpline roadSpline;

    static readonly float bridgeSegLength_Min = 15;
    static readonly float bridgeSegLength_Max = 50;
    static float ramp_midOfOpen = 0.72f; // points the the mid of the open space

    public static void OnUpdateProperties(
        float roadPartLength,
        float roadPartWidth,
        float bridgeAscentLength,
        float bridgeDescentLength,
        float crossingExtraSize,
        float rampLength,
        RoadSpline roadSpline)
    {
        TopologyHandler.roadPartLength = roadPartLength;
        TopologyHandler.roadPartWidth = roadPartWidth;
        TopologyHandler.bridgeAscentLength = bridgeAscentLength;
        TopologyHandler.bridgeDescentLength = bridgeDescentLength;
        TopologyHandler.crossingExtraSize = crossingExtraSize;
        TopologyHandler.rampLength = rampLength;
        TopologyHandler.roadSpline = roadSpline;
    }

    static void Initialize()
    {
        intersectingPercents = CalculateIntersectingPercents();
        fourPercCrossings = new();
    }

    public static void SetTopologies(List<TopologyRange> topologies, List<float> fourPercCrossings_1D)
    {
        Initialize();

        TopologyHandler.topologies = new(topologies);
        if (fourPercCrossings_1D != null)
        {
            fourPercCrossings = new();
            for (int i = 0; i < fourPercCrossings_1D.Count / 4; i++)
            {
                fourPercCrossings.Add(new()
                {
                    fourPercCrossings_1D[i * 4 + 0],
                    fourPercCrossings_1D[i * 4 + 1],
                    fourPercCrossings_1D[i * 4 + 2],
                    fourPercCrossings_1D[i * 4 + 3]
                });
            }
        }
    }

    public static void GenerateTopologies()
    {
        Initialize();

        topologies = EmptyTopologyList();

        //Debug.Log(intersectingPercents.ExtendedToString());

        for (int index = 0; index < intersectingPercents.Count; index++)
        {
            var intersection = intersectingPercents[index];
            float a = intersection.start;
            float b = intersection.end;

            bool a_crossingPossible = TopologyPossible(a, () => Crossing_GetStartEnd(a, b));
            bool b_crossingPossible = TopologyPossible(b, () => Crossing_GetStartEnd(b, a));
            bool crossing_possible = a_crossingPossible && b_crossingPossible;

            bool a_minBridge_possible = TopologyPossible(a, () => MinBridge_GetStartEnd(a));
            bool b_minBridge_possible = TopologyPossible(b, () => MinBridge_GetStartEnd(b));

            bool a_ramp_possible = TopologyPossible(a, () => Ramp_GetStartEnd(a));
            bool b_ramp_possible = TopologyPossible(b, () => Ramp_GetStartEnd(b));

            //Debug.Log(AnotherIntersectionInRange(b_bridgeStart, b_bridgeEnd, index));
            //Debug.Log("ramp possible == " + a_ramp_possible + " around " + a);
            //Debug.Log("ramp possible == " + b_ramp_possible + " around " + b);
            //Debug.Log(Ramp_GetStartEnd(a));
            //Debug.Log("Prev to " + a + ": " + GetPrevIntersecPerc(a));
            //Debug.Log("Next to " + a + ": " + GetNextIntersecPerc(a));

            var topsAndPossibles = new List<Tuple<string, bool>>()
            {
                new("crossing", crossing_possible),
                new("crossing", crossing_possible),
                new("crossing", crossing_possible),
                new("bridge_a", a_minBridge_possible),
                new("bridge_b", b_minBridge_possible),
                new("ramp_a", a_ramp_possible),
                new("ramp_b", b_ramp_possible)
            };

            List<string> topsAndPossibles_possibles = topsAndPossibles
                .Where(tuple => tuple.Item2)
                .Select(tuple => tuple.Item1)
                .ToList();
            if (topsAndPossibles_possibles.Count < 1)
            {
                Debug.LogWarning("Nothing possible -- Adding a crossing anyways.");
                AddCrossingAround(a, b, topologies);
            }


            var randomIndex = UnityEngine.Random.Range(0, topsAndPossibles_possibles.Count);
            string top = topsAndPossibles_possibles[randomIndex];


            if (top == "crossing")
            {
                AddCrossingAround(a, b, topologies);
            }
            else if (top == "bridge_a")
            {
                AddBridgeAround(a, randomBridgeSegLength(a), topologies);
            }
            else if (top == "bridge_b")
            {
                AddBridgeAround(b, randomBridgeSegLength(b), topologies);
            }
            else if (top == "ramp_a")
            {
                AddRampAround(a, topologies);
            }
            else if (top == "ramp_b")
            {
                AddRampAround(b, topologies);
            }
            else
            {
                Debug.LogWarning("Something went wrong!");
            }

            LogTopologies();
        }

        float DstTo_prev(float perc)
        {
            var prev = GetPrevIntersecPerc(perc);
            return percDistAtoNextB(prev, perc);
        }

        float DstTo_next(float perc)
        {
            var next = GetNextIntersecPerc(perc);
            return percDistAtoNextB(perc, next);
        }

        bool TopologyPossible(float perc, Func<(float, float)> GetStartEnd)
        {
            var dstTo_prev = DstTo_prev(perc);
            var dstTo_next = DstTo_next(perc);

            var (top_start, top_end) = GetStartEnd();
            var (dstTo_top_start, dstTo_top_end) = (percDistAtoNextB(top_start, perc), percDistAtoNextB(perc, top_end));
            return (dstTo_top_start < dstTo_prev) && (dstTo_top_end < dstTo_next);
        }

        float randomBridgeSegLength(float perc)
        {
            float dstTo_prev = DstTo_prev(perc);
            float dstTo_next = DstTo_next(perc);

            float maxPercDst = dstTo_next + dstTo_prev;

            float extremititesPercDst = (bridgeAscentLength + bridgeDescentLength) / splineLength;
            float minBridgePercDst = extremititesPercDst + bridgeSegLength_Min / splineLength;
            float maxBridgePercDst = extremititesPercDst + bridgeSegLength_Max / splineLength;
            maxBridgePercDst = Mathf.Min(maxPercDst, maxBridgePercDst);

            float bridgePercLength = UnityEngine.Random.Range(minBridgePercDst, maxBridgePercDst);
            float bridgeSegLength = (bridgePercLength - extremititesPercDst) * splineLength;

            if (bridgeSegLength < 0)
            {
                Debug.LogWarning("bridgeSegLength < 0");
            }

            return bridgeSegLength;
        }

        float percDistAtoNextB(float a, float b) // a must be before b! This needs to be specified, as the spline loops and has no sense of before / after on its own
        {
            float dst = b - a;
            if (dst < 0) dst += 1;
            return dst;
        }

        (float start, float end) MinBridge_GetStartEnd(float perc)
        {
            float bridgeSegLength_Min_Perc = bridgeSegLength_Min / splineLength;
            float bridgeAscentLength_Perc = bridgeAscentLength / splineLength;
            float bridgeDescentLength_Perc = bridgeDescentLength / splineLength;

            float start = perc - (bridgeSegLength_Min_Perc * 0.5f + bridgeAscentLength_Perc);
            float end = perc + (bridgeSegLength_Min_Perc * 0.5f + bridgeDescentLength_Perc);

            return StartEnd01(start, end);
        }

        float GetNextIntersecPerc(float perc)
        {
            float next = intersectingPercents_Flat().OrderBy(otherPerc =>
            {
                if (otherPerc <= perc) return otherPerc + 1; // <= instead of < to make own perc the farthest
                return otherPerc;
            }).First();

            return next;
        }

        float GetPrevIntersecPerc(float perc)
        {
            float prev = intersectingPercents_Flat().OrderBy(otherPerc =>
            {
                if (otherPerc >= perc) return otherPerc - 1; // >= instead of > to make own perc the farthest
                return otherPerc;
            }).Last();

            return prev;
        }

        List<float> intersectingPercents_Flat()
        {
            List<float> allIntersections = new();
            for (int i = 0; i < intersectingPercents.Count; i++)
            {
                var intersection = intersectingPercents[i];
                allIntersections.AddRange(new List<float>() { intersection.Item1, intersection.Item2 });
            }
            return allIntersections;
        }

        //topologies.ForEach(x => Debug.Log(x.Item1 + ", " + x.Item2));

    }

    public static void OnSplineChanged()
    {
        int intersectionCountBefore = intersectingPercents.Count;

        intersectingPercents = CalculateIntersectingPercents();

        var newTopologies = new List<Tuple<float, TopologyType>>();

        if (intersectionCountBefore != intersectingPercents.Count)
        {
            Debug.LogWarning("Not Implemented yet!");
            return;
        }

        Debug.LogWarning("Not Implemented yet!");
    }

    static int GetTopologyIndex(float t, List<TopologyRange> topologies)
    {
        for (int i = 0; i < topologies.Count; i++)
        {
            float topologyEndT = topologies[i].percEnd;
            if (t < topologyEndT)
            {
                return i;
            }
        }
        return topologies.Count > 1 ? (topologies.Count - 1) : 0;
    }

    static (TopologyType type, float tInSegment) GetTopologyAndTInTopology(float dstAlongSpline)
    {
        return GetTopologyAndTInTopology(dstAlongSpline, topologies);
    }

    public static (TopologyType type, float tInSegment) GetTopologyAndTInTopology(float dstAlongSpline, List<TopologyRange> topologies)
    {
        return GetTopologyAndTInTopology(dstAlongSpline, splineLength, topologies);
    }

    public static (TopologyType type, float tInSegment) GetTopologyAndTInTopology(float dstAlongSpline, float splineLength, List<TopologyRange> topologies)
    {
        float perc = dstAlongSpline / splineLength;

        int index = GetTopologyIndex(perc, topologies);
        float topologyEndT = topologies[index].percEnd;
        float topologyStartT = index > 0 ? topologies[index - 1].percEnd : 0;
        float tInTopology = (perc - topologyStartT) / (topologyEndT - topologyStartT);
        return (topologies[index].type, tInTopology);
    }



    static void AddBridgeAround(float perc, float midSectionLength, List<TopologyRange> topologies)
    {
        float tRange_bridge = midSectionLength / splineLength;
        float tRange_ascent = bridgeAscentLength / splineLength;
        float tRange_descent = bridgeDescentLength / splineLength;

        float t_ascentEnd = perc - tRange_bridge * 0.5f;
        float t_ascentStart = t_ascentEnd - tRange_ascent;
        float t_bridgeEnd = perc + tRange_bridge * 0.5f;
        float t_descentEnd = t_bridgeEnd + tRange_descent;

        (t_ascentStart, t_ascentEnd) = StartEnd01(t_ascentStart, t_ascentEnd);
        (t_ascentEnd, t_bridgeEnd) = StartEnd01(t_ascentEnd, t_bridgeEnd);
        (t_bridgeEnd, t_descentEnd) = StartEnd01(t_bridgeEnd, t_descentEnd);

        AddTopologyToList(t_ascentStart, t_ascentEnd, TopologyType.BridgeAscent, ref topologies);
        AddTopologyToList(t_ascentEnd, t_bridgeEnd, TopologyType.Bridge, ref topologies);
        AddTopologyToList(t_bridgeEnd, t_descentEnd, TopologyType.BridgeDescent, ref topologies);
    }

    static void AddCrossingAround(float a_perc, float b_perc, List<TopologyRange> topologies)
    {
        var (a_start, a_end) = Crossing_GetStartEnd(a_perc, b_perc);
        var (b_start, b_end) = Crossing_GetStartEnd(b_perc, a_perc);

        AddTopologyToList(a_start, a_end, TopologyType.Crossing, ref topologies);
        AddTopologyToList(b_start, b_end, TopologyType.Crossing, ref topologies);

        List<float> newCrossing;
        var (_, a_segIndex, a_t) = roadSpline.CalculatePointAtDist(a_perc * splineLength);
        Vector2 a_tangent = roadSpline.UnitTangent(a_segIndex, a_t);
        var (_, b_segIndex, b_t) = roadSpline.CalculatePointAtDist(b_perc * splineLength);
        Vector2 b_tangent = roadSpline.UnitTangent(b_segIndex, b_t);
        if (Vector2.SignedAngle(a_tangent, b_tangent) > 0)
            newCrossing = new() { a_start, b_end, a_end, b_start };
        else
            newCrossing = new() { a_start, b_start, a_end, b_end };
        fourPercCrossings.Add(newCrossing);
    }

    static (float start, float end) Crossing_GetStartEnd(float a_perc, float b_perc)
    {
        float splineLength = roadSpline.TotalLength;
        (Vector2 _, int segIndex, float t) = roadSpline.CalculatePointAtDist(a_perc * splineLength);
        (Vector2 _, int segIndex_other, float t_other) = roadSpline.CalculatePointAtDist(b_perc * splineLength);
        Vector2 a_tangent = roadSpline.UnitTangent(segIndex, t);
        Vector2 b_tangent = roadSpline.UnitTangent(segIndex_other, t_other);

        float angle = Vector2.Angle(a_tangent, b_tangent) * Mathf.Deg2Rad;    // [0, Mathf.PI)
        if (angle > Mathf.PI * 0.5f)
            angle = Mathf.PI - angle;                                           // [0, Mathf.PI * 0.5f)
        float invAngle = Mathf.PI * 0.5f - angle;                                      // Inversed
        float tan = Mathf.Tan(invAngle);
        float midToEdge = 0.5f * roadPartWidth * (1 + tan);
        float unitRange = 2 * (midToEdge + crossingExtraSize * roadPartWidth);
        float percRange = unitRange / splineLength;
        float percStart = a_perc - percRange * 0.5f;
        float percEnd = a_perc + percRange * 0.5f;

        //Debug.Log("Crossing start end around " + a_perc + " and " + b_perc + ": (" + tStart + ", " + tEnd + ")");
        return StartEnd01(percStart, percEnd);
    }

    /// <summary>
    /// Returns an estimate of the start-end-perc-range according to Crossing_GetStartEnd() where neither the tangents nor the road width is known
    /// </summary>
    /// <returns></returns>
    public static float Crossing_EstimatedStartEndRange()
    {
        const float estimatedRoadPartWidth = 5;
        float splineLength = roadSpline.TotalLength;

        float midToEdge = 0.5f * estimatedRoadPartWidth;
        float unityRange = 2 * (midToEdge + crossingExtraSize * roadPartWidth);
        float percRange = unityRange / splineLength;
        return percRange;
    }

    static void AddRampAround(float perc, List<TopologyRange> topologies)
    {
        var (perc_start, perc_end) = Ramp_GetStartEnd(perc);
        AddTopologyToList(perc_start, perc_end, TopologyType.Ramp, ref topologies);
    }

    static (float start, float end) Ramp_GetStartEnd(float perc)
    {
        float ramp_perc = rampLength / splineLength;
        float start = perc - (ramp_perc * ramp_midOfOpen);
        float end = perc + (ramp_perc * (1 - ramp_midOfOpen));
        return StartEnd01(start, end);
    }

    public static (float start, float end) StartEnd01(float start, float end)
    {
        if (start < 0) start += 1;
        if (end > 1) end -= 1;
        return (start, end);
    }

    /// <param name="tStart">MUST BE 01</param>
    /// <param name="tEnd">MUST BE 01</param>
    public static void AddTopologyToList(float tStart, float tEnd, TopologyType topType, ref List<TopologyRange> topologies)
    {
        if (!(tEnd >= 0 && tEnd <= 1) ||
            !(tStart >= 0 && tStart <= 1))
        {
            Debug.LogWarning($"ERROR: Range is not 01-normed: ({tStart}, {tEnd})");
            return;
        }
        if (tEnd - tStart > 1) // To catch recursion
        {
            Debug.LogWarning($"ERROR: Provided an invalid tRange: ({tStart}, {tEnd})");
            return;
        }

        //Debug.Log($"Adding Topology with ({tStart}, {tEnd})");

        if (tEnd > tStart) // Usual case, no wrapping around
        {
            for (int i = 0; i < topologies.Count; i++)
            {
                if (topologies[i].percEnd >= tEnd)
                {
                    bool noTopologyBetween = (i > 1 && topologies[i - 1].percEnd == tStart) || (i == 0 && tStart == 0);
                    bool needToSplitTop = !noTopologyBetween;
                    if (needToSplitTop)
                    {
                        topologies.InsertRange(i, new List<TopologyRange>() {
                        new(tStart, topologies[i].type),
                        new(tEnd, topType)
                    });
                    }
                    else
                    {
                        topologies.Insert(i, new(tEnd, topType));
                    }

                    if (tEnd == 1)
                    {
                        topologies.RemoveAt(topologies.Count - 1);
                    }

                    break;
                }
            }
        }
        else // Wrapping around
        {
            //LogTopologies();
            for (int i = 0; i < topologies.Count; i++)
            {
                if (topologies[i].percEnd > tEnd)
                {
                    topologies.RemoveRange(0, i);
                    topologies.Insert(0, new(tEnd, topType));
                    break;
                }
            }
            //LogTopologies();

            for (int i = 0; i < topologies.Count; i++)
            {
                if (topologies[i].percEnd > tStart)
                {
                    topologies.RemoveRange(i, topologies.Count - i);// Alle weiteren entfernen

                    topologies.Add(new(tStart, TopologyType.Normal)); // Letzte Topology darf nicht verlängert werden -> Normal einfügen

                    topologies.Add(new(1, topType)); // Topology einfügen

                    break;
                }
            }
        }

        //LogTopologies();
    }

    public static List<(float a, float b)> CalculateIntersectingPercents()
    {
        List<(float a, float b)> intersectingPercents = new();

        var intersections = roadSpline.CalculateIntersections();
        //Debug.Log(intersections.ExtendedToString());

        splineLength = roadSpline.TotalLength;

        foreach (var intersection in intersections)
        {
            float percentOne = roadSpline.GetDstAlongSpline(intersection.Item1.Item1, intersection.Item1.Item2) / splineLength;
            float percentTwo = roadSpline.GetDstAlongSpline(intersection.Item2.Item1, intersection.Item2.Item2) / splineLength;
            //Debug.Log("percentOne: " + percentOne + ", percentTwo: " + percentTwo);

            intersectingPercents.Add(new(percentOne, percentTwo));

            // All positions must be equal
            //Place(roadSpline.CalculatePointAtT(intersection.Item1.Item1, intersection.Item1.Item2));
            //Place(roadSpline.CalculatePointAfterRef(percentOne * splineLength, 0, 0, 10f).Item1);
            //Place(roadSpline.CalculatePointAtT(intersection.Item2.Item1, intersection.Item2.Item2));
            //Place(roadSpline.CalculatePointAfterRef(percentTwo * splineLength, 0, 0, 10f).Item1);
            //void Place(Vector2 pt)
            //{
            //    GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
            //    cube.transform.position = new(pt.x, 0, pt.y);
            //}
        }

        return intersectingPercents;
    }

    public static List<TopologyRange> EmptyTopologyList()
    {
        return new()
        {
            new(1, TopologyType.Normal)
        };
    }

    static void LogTopologies()
    {
        LogTopologies(topologies);
    }

    public static void LogTopologies(List<TopologyRange> topologies)
    {
        if (topologies == null)
        {
            Debug.LogWarning("Please first generate topologies.");
            return;
        }

        string output = "";
        topologies.ForEach(x => output += x + "\n");
        Debug.Log(output);
    }

    static float LoopT(float t)
    {
        return t > 1 ? t - 1 : t < 0 ? t + 1 : t;
    }
}

[Serializable]
public class TopologyRange
{
    [SerializeField]
    public float percEnd;
    [SerializeField]
    public TopologyType type;

    public TopologyRange(float percEnd, TopologyType type)
    {
        this.percEnd = percEnd;
        this.type = type;
    }

    public override string ToString()
    {
        return "endT: " + percEnd + ", type: " + type;
    }
}

[Serializable]
public enum TopologyType
{
    Normal,
    BridgeAscent,
    Bridge,
    BridgeDescent,
    Ramp,
    Crossing
}