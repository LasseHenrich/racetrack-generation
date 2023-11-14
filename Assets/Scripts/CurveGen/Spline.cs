using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class Spline
{
    public ObservableList<Vector2> controls;
    protected bool isClosed;
    bool autoSetControlPoints;
    List<float> _segmentLengths;

    public bool IsClosed { get { return isClosed; } }

    public List<Vector2[]> pointsInSegment; // Apparently, we call GetPointsInSegment() a LOT of times, so the garbage collector struggled a bit. To not let the garbage be collected, we created this variable here.

    public float GetSegmentLength(int segIndex)
    {
        if (_segmentLengths == null || segIndex >= _segmentLengths.Count)
            CalculateAllSegmentLengths();
        return _segmentLengths[segIndex];
    }


    public Spline(Vector2 centre)
    {
        controls = new(this);
        controls.AddRange(new List<Vector2>()
        {
            centre + Vector2.left * 10,
            centre + (Vector2.left + Vector2.up) * 5f,
            centre + (Vector2.right + Vector2.down) * 5f,
            centre + Vector2.right * 10
        });
    }

    public Spline(List<Vector2> newCtrls, bool controlPointsIncluded = false, bool alreadyClosed = false)
    {
        if (!controlPointsIncluded)
        {
            this.controls = new(this);
            for (int i = 0; i < newCtrls.Count; i++)
            {
                if (i > 0)
                    this.controls.Add(newCtrls[i] - Vector2.one);
                this.controls.Add(newCtrls[i]);
                if (i < newCtrls.Count - 1)
                    this.controls.Add(newCtrls[i] + Vector2.one);
            }
        }
        else
        {
            this.controls = new(this);
            this.controls.AddRange(newCtrls);

            if (alreadyClosed)
                isClosed = true;
        }
    }

    #region Editor

    public bool AutoSetControlPoints
    {
        get
        {
            return autoSetControlPoints;
        }
        set
        {
            if (autoSetControlPoints != value)
            {
                autoSetControlPoints = value;
                if (autoSetControlPoints)
                {
                    AutoSetAllControlPoints();
                }
            }
        }
    }

    public Vector2 this[int i]
    {
        get
        {
            return controls[i];
        }
    }

    public int NumPoints
    {
        get
        {
            return controls.Count;
        }
    }

    public int NumSegments
    {
        get
        {
            return controls.Count / 3; // same as [(points.Count - 4) / 3 + 1 + ((isClosed)?1:0)] rounded down
        }
    }

    public int SegmentOfControl(int control)
    {
        return control / 3;
    }

    public void AddSegment(Vector2 anchorPos)
    {
        controls.Add(controls[controls.Count - 1] * 2 - controls[controls.Count - 2]);
        controls.Add((controls[controls.Count - 1] + anchorPos) * .5f);
        controls.Add(anchorPos);

        if (autoSetControlPoints)
        {
            AutoSetAllAffectedControlPoints(controls.Count - 1);
        }
    }

    public virtual void SplitSegment(Vector2 anchorPos, int segmentIndex)
    {
        controls.InsertRange(segmentIndex * 3 + 2, new Vector2[] { Vector2.zero, anchorPos, Vector2.zero });
        if (autoSetControlPoints)
        {
            AutoSetAllAffectedControlPoints(segmentIndex * 3 + 3);
        }
        else
        {
            AutoSetAnchorControlPoints(segmentIndex * 3 + 3);
        }
    }


    public void DeleteSegment(int anchorIndex)
    {
        if (NumSegments > 2 || !isClosed && NumSegments > 1)
        {
            if (anchorIndex == 0)
            {
                if (isClosed)
                {
                    controls[controls.Count - 1] = controls[2];
                }
                controls.RemoveRange(0, 3);
            }
            else if (anchorIndex == controls.Count - 1 && !isClosed)
            {
                controls.RemoveRange(anchorIndex - 2, 3);
            }
            else
            {
                controls.RemoveRange(anchorIndex - 1, 3);
            }
        }
    }

    public void CalculateAllPointsInSegments()
    {
        pointsInSegment = new();
        for (int i = 0; i < NumSegments; i++)
        {
            pointsInSegment.Add(new Vector2[] { controls[i * 3], controls[i * 3 + 1], controls[i * 3 + 2], controls[LoopIndex(i * 3 + 3)] });
        }
    }

    public Vector2[] GetPointsInSegment(int segIndex)
    {
        return pointsInSegment[segIndex];
    }

    public void MovePoint(int i, Vector2 pos)
    {
        Vector2 deltaMove = pos - controls[i];
        controls[i] = pos;

        if (i % 3 == 0)
        {
            if (i + 1 < controls.Count || isClosed)
            {
                controls[LoopIndex(i + 1)] += deltaMove;
            }
            if (i - 1 >= 0 || isClosed)
            {
                controls[LoopIndex(i - 1)] += deltaMove;
            }
        }
        else
        {
            bool nextPointIsAnchor = (i + 1) % 3 == 0;
            int correspondingControlIndex = (nextPointIsAnchor) ? i + 2 : i - 2;
            int anchorIndex = (nextPointIsAnchor) ? i + 1 : i - 1;

            if (correspondingControlIndex >= 0 && correspondingControlIndex < controls.Count || isClosed)
            {
                float dst = (controls[LoopIndex(anchorIndex)] - controls[LoopIndex(correspondingControlIndex)]).magnitude;
                Vector2 dir = (controls[LoopIndex(anchorIndex)] - pos).normalized;
                controls[LoopIndex(correspondingControlIndex)] = controls[LoopIndex(anchorIndex)] + dir * dst;
            }
        }
    }

    public void AutoSetAllAffectedControlPoints(int updatedAnchorIndex) // AutoSet the control Points of both affected anchors
    {
        for (int i = updatedAnchorIndex - 3; i <= updatedAnchorIndex; i += 3)
        {
            if (i >= 0 && i < controls.Count)
            {
                AutoSetAnchorControlPoints(i);
            }
        }
    }

    public void AutoSetAllControlPoints()
    {
        for (int i = 0; i < controls.Count; i += 3)
        {
            AutoSetAnchorControlPoints(i);
        }
    }

    void AutoSetAnchorControlPoints(int anchorIndex)
    {
        Vector2 anchorPos = controls[anchorIndex];
        //Debug.Log(anchorIndex);
        Vector2 dir = Vector2.zero;
        float[] neighbourDistances = new float[2];

        if (anchorIndex - 3 >= 0 || isClosed)
        {
            Vector2 offset = controls[LoopIndex(anchorIndex - 3)] - anchorPos;
            dir += offset.normalized;
            neighbourDistances[0] = offset.magnitude;
            //Debug.Log(offset);
        }
        if (anchorIndex + 3 >= 0 || isClosed)
        {
            Vector2 offset = controls[LoopIndex(anchorIndex + 3)] - anchorPos;
            dir -= offset.normalized;
            neighbourDistances[1] = -offset.magnitude;
            //Debug.Log(offset);
        }

        dir.Normalize();

        for (int i = 0; i < 2; i++)
        {
            int controlIndex = anchorIndex + i * 2 - 1; // anchorIndex - 1 || anchorIndex + 1
            if (controlIndex >= 0 && controlIndex < controls.Count || isClosed)
            {
                controls[LoopIndex(controlIndex)] = anchorPos + dir * neighbourDistances[i] * 0.5f;
                //Mess.instance.InstantiateAtPos(points[controlIndex], controlIndex.ToString());
            }
        }
        //Debug.Log("anchorPos" + anchorIndex + ": " + anchorPos);
    }

    public void ToggleClosed()
    {
        isClosed = !isClosed;

        if (isClosed)
        {
            controls.Add(controls[controls.Count - 1] * 2 - controls[controls.Count - 2]);
            controls.Add(controls[0] * 2 - controls[2]);
            if (autoSetControlPoints)
            {
                AutoSetAnchorControlPoints(0);
                AutoSetAnchorControlPoints(controls.Count - 3);
            }
        }
        else
        {
            controls.RemoveRange(controls.Count - 2, 2);
        }
    }

    public int LoopIndex(int i) => (i + controls.Count) % controls.Count;

    public int LoopIndex<T>(int i, List<T> points) => (i + points.Count) % points.Count;

    public int LoopIndex(int i, int count) => (i + count) % count;


    #endregion

    #region Runtime Calculations

    /// <summary>
    /// Evenly spaced points on whole path
    /// </summary>
    /// <param name="spacing"></param>
    /// <param name="resolution"></param>
    public List<Vector2> CalculateEvenlySpacedPoints(float spacing, float resolution = 1)
    {
        return CalculateEvenlySpacedPoints(spacing, out _, resolution);
    }

    public List<Vector2> CalculateEvenlySpacedPoints(float spacing, out List<int> indicesOfSegments, float resolution = 1)
    {
        return CalculateEvenlySpacedPoints(spacing, out indicesOfSegments, out _, resolution);
    }

    /// <summary>
    /// Evenly spaced points on whole path
    /// </summary>
    /// <param name="spacing"></param>
    /// <param name="indicesOfSegments">Points to the index of the first point in the segment</param>
    /// <param name="resolution"></param>
    public List<Vector2> CalculateEvenlySpacedPoints(float spacing, out List<int> indicesOfSegments, out List<float> ts, float resolution = 1)
    {
        float totalLength = 0;
        for (int segmentIndex = 0; segmentIndex < NumSegments; segmentIndex++)
            totalLength += GetSegmentLength(segmentIndex);
        int pointCount = (int)(totalLength / spacing);
        spacing = totalLength / pointCount;

        List<Vector2> evenlySpacedPoints = new List<Vector2>();
        Vector2 previousPoint = controls[0];
        float dstSinceLastEvenPoint = 0; // distance since last evenly spaced point

        indicesOfSegments = new();
        ts = new();

        for (int segmentIndex = 0; segmentIndex < NumSegments; segmentIndex++)
        {
            indicesOfSegments.Add(evenlySpacedPoints.Count);

            Vector2[] p = GetPointsInSegment(segmentIndex);
            float estimatedSegLength = GetSegmentLength(segmentIndex);
            int divisions = Mathf.CeilToInt(estimatedSegLength * resolution * 10);
            float t = 0f;
            while (t <= 1f)
            {
                Vector2 pointOnCurve = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], t);
                dstSinceLastEvenPoint += Vector2.Distance(previousPoint, pointOnCurve);

                while (dstSinceLastEvenPoint >= spacing)
                {
                    float overshootDst = dstSinceLastEvenPoint - spacing;
                    Vector2 newEvenlySpacedPoint = pointOnCurve + (previousPoint - pointOnCurve).normalized * overshootDst;
                    evenlySpacedPoints.Add(newEvenlySpacedPoint);
                    dstSinceLastEvenPoint = overshootDst;
                    previousPoint = newEvenlySpacedPoint;

                    ts.Add(t - (overshootDst / estimatedSegLength));
                }

                previousPoint = pointOnCurve;
                t += 1f / divisions;
            }
        }

        if (dstSinceLastEvenPoint < spacing * 0.99f)
        {
            evenlySpacedPoints.RemoveAt(evenlySpacedPoints.Count - 1);
            float maxMove = dstSinceLastEvenPoint;
            for (int i = 1; i < evenlySpacedPoints.Count; i++)
            {
                int nextPoint = (i + 1) % evenlySpacedPoints.Count;
                Vector2 direction = (evenlySpacedPoints[nextPoint] - evenlySpacedPoints[i]).normalized;
                Vector2 move = direction * maxMove * ((float)i / (evenlySpacedPoints.Count - 1));
                evenlySpacedPoints[i] += move;
            }
        }

        return evenlySpacedPoints;
    }

    /// <summary>
    /// Evenly spaced point on one segment
    /// </summary>
    /// <param name="segIndex"></param>
    /// <param name="spacing"></param>
    /// <param name="previousPoint"></param>
    /// <param name="resolution"></param>
    /// <returns></returns>
    public List<Vector2> CalculateEvenlySpacedPointsInSegment(int segIndex, float spacing, out List<float> ts, float resolution = 1)
    {
        Vector2[] p = GetPointsInSegment(segIndex);

        Vector2 previousPoint = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], 0);

        float estimatedSegLength = GetSegmentLength(segIndex);
        int divisions = Mathf.CeilToInt(estimatedSegLength * resolution * 10);
        float dstSinceLastEvenPoint = 0; // distance since last evenly spaced point

        List<Vector2> evenlySpacedPoints = new();
        ts = new();

        float t = 0f;

        while (t <= 1f)
        {
            t += 1f / divisions;
            Vector2 pointOnCurve = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], t);
            dstSinceLastEvenPoint += Vector2.Distance(previousPoint, pointOnCurve);

            while (dstSinceLastEvenPoint >= spacing)
            {
                float overshootDst = dstSinceLastEvenPoint - spacing;
                Vector2 newEvenlySpacedPoint = pointOnCurve + (previousPoint - pointOnCurve).normalized * overshootDst;
                evenlySpacedPoints.Add(newEvenlySpacedPoint);
                dstSinceLastEvenPoint = overshootDst;
                previousPoint = newEvenlySpacedPoint;

                ts.Add(t - (overshootDst / estimatedSegLength));
            }

            previousPoint = pointOnCurve;
        }

        return evenlySpacedPoints;
    }

    public void CalculateAllSegmentLengths()
    {
        //Debug.Log("Calculating all segment lengths");
        _segmentLengths = new();
        for (int i = 0; i < NumSegments; i++)
            CalculateSegmentLength(i);
    }

    // [Primer, §25]
    void CalculateSegmentLength(int segIndex)
    {
        // This does seem to work in MOST scenarios, but not all
        //Vector2[] p = GetPointsInSegment(segIndex);
        //float controlNetLength = Vector2.Distance(p[0], p[1]) + Vector2.Distance(p[1], p[2]) + Vector2.Distance(p[2], p[3]);
        //float estimatedCurveLength = Vector2.Distance(p[0], p[3]) + controlNetLength / 2f;
        //return estimatedCurveLength;

        int divisions = 50;

        Vector2[] p = GetPointsInSegment(segIndex);

        Vector2 previousPoint = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], 0);

        float length = 0;

        for (int i = 0; i < divisions; i++)
        {
            float t = (i + 1) / (float)divisions;
            var point = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], t);
            length += Vector2.Distance(point, previousPoint);
            previousPoint = point;
        }

        if (segIndex < _segmentLengths.Count)
            _segmentLengths[segIndex] = length;
        else if (segIndex == _segmentLengths.Count)
            _segmentLengths.Add(length);
        else
        {
            Debug.LogWarning("ERROR: Something went wrong!");
        }
    }

    public float TotalLength
    {
        get
        {
            float totalLength = 0;
            for (int segmentIndex = 0; segmentIndex < NumSegments; segmentIndex++)
                totalLength += GetSegmentLength(segmentIndex);
            return totalLength;
        }
    }

    public Vector2[] CalculateForwards(int index, Vector2[] points, float[] t)
    {
        Vector2[] forwards = new Vector2[t.Length];
        Vector2[] pThis = GetPointsInSegment(index);
        Vector2[] pLast;
        Vector2[] pNext;
        float gap = 0.0001f;

        float thisSegmentLength;
        float lastSegmentLength;
        float nextSegmentLength;

        pLast = GetPointsInSegment(LoopIndex(index - 1, NumSegments));
        thisSegmentLength = GetSegmentLength(index);
        lastSegmentLength = GetSegmentLength(LoopIndex(index - 1, NumSegments));

        pNext = GetPointsInSegment(LoopIndex(index + 1, NumSegments));
        nextSegmentLength = GetSegmentLength(LoopIndex(index + 1, NumSegments));

        for (int i = 0; i < t.Length; i++)
        {
            // With a previous point
            if (t[i] - gap < 0) // Last curve
                forwards[i] += points[i] - xy0ToXY(Bezier.EvaluateCubic(pLast[0], pLast[1], pLast[2], pLast[3], ((t[i] - gap) * (thisSegmentLength / lastSegmentLength)) + 1)); // t[i] is relative: has to be recalculated to fit last segment
            else if (t[i] - gap <= 1) // This curve
                forwards[i] += points[i] - xy0ToXY(Bezier.EvaluateCubic(pThis[0], pThis[1], pThis[2], pThis[3], t[i] - gap));
            else // Next curve
                forwards[i] += points[i] - xy0ToXY(Bezier.EvaluateCubic(pNext[0], pNext[1], pNext[2], pNext[3], (((t[i] - 1) - gap) * (thisSegmentLength / nextSegmentLength))));

            // With a next point
            if (t[i] + gap > 1) // Next curve
                forwards[i] += xy0ToXY(Bezier.EvaluateCubic(pNext[0], pNext[1], pNext[2], pNext[3], (((t[i] - 1) + gap) * (thisSegmentLength / nextSegmentLength)))) - points[i];
            else if (t[i] + gap >= 0) // This curve
                forwards[i] += xy0ToXY(Bezier.EvaluateCubic(pThis[0], pThis[1], pThis[2], pThis[3], t[i] + gap)) - points[i];
            else // Last curve
                forwards[i] += xy0ToXY(Bezier.EvaluateCubic(pLast[0], pLast[1], pLast[2], pLast[3], ((t[i] + gap) * (thisSegmentLength / lastSegmentLength)) + 1)) - points[i];

            forwards[i].Normalize();
        }

        Vector2 xy0ToXY(Vector3 xy0)
        {
            return new Vector2(xy0.x, xy0.y);
        }

        return forwards;
    }

    public Vector2 CalculatePointAtT(int segIndex, float t)
    {
        Vector2[] p = GetPointsInSegment(segIndex);
        var p3D = Bezier.EvaluateCubic(p[0], p[1], p[2], p[3], t);
        return new(p3D.x, p3D.y);
    }

    //public float closestTonSegment(int segIndex, Vector2 point)
    //{
    //    var evenlySpacedPoints = CalculateEvenlySpacedPointsInSegment(segIndex, 0.1f);
    //
    //}

    // Invserse Function of CalculatePointAtDist
    // Note that overshooting should also be considered (but doesn't really matter as we're taking very small steps)
    public float GetDstAlongSpline(int segIndex, float t, bool highResolution = true)
    {
        float length = 0;
        for (int i = 0; i < segIndex; i++)
        {
            length += GetSegmentLength(i);
        }

        CalculateEvenlySpacedPointsInSegment(segIndex, highResolution ? 0.01f : 0.1f, out var ts);
        int t_index = 0;
        for (int i = 0; i < ts.Count; i++)
        {
            if (t < ts[i])
            {
                t_index = i;
                break;
            }
        }
        float percentageAlongSeg = (float)t_index / ts.Count;
        float lengthInSegment = percentageAlongSeg * GetSegmentLength(segIndex);
        length += lengthInSegment;

        return length;
    }

    #endregion

    #region Bounding Box

    public (Vector2, Vector2) BoundingBox(int segmentIndex)
    {
        List<Vector2> points = GetPointsInSegment(segmentIndex).ToList();
        return BoundingBox(points);
    }

    public (Vector2, Vector2) BoundingBox(List<Vector2> ps)
    {
        var p0 = ps[0]; var p1 = ps[1]; var p2 = ps[2]; var p3 = ps[3];

        List<float> roots = FindRoots(p0, p1, p2, p3);

        for (int i = 0; i < roots.Count; i++)
        {
            if (roots[i] < 0 || roots[i] > 1)
                roots.RemoveAt(i);
        }
        roots.AddRange(new List<float>() { 0, 1 });

        List<Vector2> values = roots.Select(t => Bezier.EvaluateCubic(p0, p1, p2, p3, t)).Select(x => new Vector2(x.x, x.y)).ToList();

        float lowestX = values.Min(v => v.x);
        float lowestY = values.Min(v => v.y);
        float highestX = values.Max(v => v.x);
        float highestY = values.Max(v => v.y);

        return (new(lowestX, lowestY), new(highestX, highestY));
    }

    // See https://snoozetime.github.io/2018/05/22/bezier-curve-bounding-box.html
    private static List<float> FindRoots(Vector2 p0, Vector2 p1, Vector2 p2, Vector2 p3)
    {
        Vector2 a = 3 * (-p0 + 3 * p1 - 3 * p2 + p3);
        Vector2 b = 6 * (p0 - 2 * p1 + p2);
        Vector2 c = 3 * (p1 - p0);

        List<float> roots = new List<float>();

        // along x
        float discriminantX = b.x * b.x - 4 * a.x * c.x;
        if (discriminantX < 0)
        {
            // No roots
        }
        else if (discriminantX == 0)
        {
            // one real root
            float rootx = (-b.x) / (2 * a.x);
            if (rootx >= 0 && rootx <= 1)
            {
                roots.Add(rootx);
            }
        }
        else if (discriminantX > 0)
        {
            // Two real roots
            float rootx1 = (-b.x + Mathf.Sqrt(discriminantX)) / (2 * a.x);
            float rootx2 = (-b.x - Mathf.Sqrt(discriminantX)) / (2 * a.x);
            if (rootx1 >= 0 && rootx1 <= 1)
            {
                roots.Add(rootx1);
            }
            if (rootx2 >= 0 && rootx2 <= 1)
            {
                roots.Add(rootx2);
            }
        }

        // along y
        float discriminantY = b.y * b.y - 4 * a.y * c.y;
        if (discriminantY < 0)
        {
            // No roots
        }
        else if (discriminantY == 0)
        {
            // one real root
            float rooty = (-b.y) / (2 * a.y);
            if (rooty >= 0 && rooty <= 1)
            {
                roots.Add(rooty);
            }
        }
        else if (discriminantY > 0)
        {
            // Two real roots
            float rooty1 = (-b.y + Mathf.Sqrt(discriminantY)) / (2 * a.y);
            float rooty2 = (-b.y - Mathf.Sqrt(discriminantY)) / (2 * a.y);
            if (rooty1 >= 0 && rooty1 <= 1)
            {
                roots.Add(rooty1);
            }
            if (rooty2 >= 0 && rooty2 <= 1)
            {
                roots.Add(rooty2);
            }
        }

        return roots;
    }

    #endregion

    public override string ToString()
    {
        string output = "";
        for (int i = 0; i < controls.Count; i++)
        {
            output += controls[i];
            if (i < controls.Count - 1) output += ", ";
        }
        return output;
    }
}

public class ObservableList<T>
{
    readonly Spline observedObject; // Can't serialize this, would cause loop
    Action UpdateSegLengths { get { return observedObject.CalculateAllSegmentLengths; } }
    Action UpdatePointsInSegs { get { return observedObject.CalculateAllPointsInSegments; } }

    [SerializeField]
    List<T> myList;

    // Indexer declaration
    public T this[int index]
    {
        get
        {
            return myList[index];
        }
        set
        {
            myList[index] = value;

            UpdatePointsInSegs();
            UpdateSegLengths();
        }
    }

    public ObservableList(Spline observedObject)
    {
        this.observedObject = observedObject;
        myList = new();
    }

    public ObservableList(Spline observedObject, IEnumerable<T> collection)
    {
        this.observedObject = observedObject;
        myList = new(collection);

        UpdatePointsInSegs();
        UpdateSegLengths();
    }

    public void Add(T a)
    {
        myList.Add(a);

        UpdatePointsInSegs();
        UpdateSegLengths();
    }

    public void AddRange(ICollection<T> coll)
    {
        myList.AddRange(coll);

        UpdatePointsInSegs();
        UpdateSegLengths();
    }

    public void InsertRange(int index, ICollection<T> coll)
    {
        myList.InsertRange(index, coll);

        UpdatePointsInSegs();
        UpdateSegLengths();
    }

    public void RemoveRange(int index, int count)
    {
        myList.RemoveRange(index, count);

        UpdatePointsInSegs();
        UpdateSegLengths();
    }

    public int Count
    {
        get
        {
            return myList.Count;
        }
    }

    public List<T> GetRange(int index, int count)
    {
        return myList.GetRange(index, count);
    }

    public List<T> GetList()
    {
        return myList;
    }
}
