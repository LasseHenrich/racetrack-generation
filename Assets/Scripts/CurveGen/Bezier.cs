using UnityEngine;

public static class Bezier
{
    private static Vector3 EvaluateQuadratic(Vector3 a, Vector3 b, Vector3 c, float t)
    {
        Vector3 p0 = Vector3.Lerp(a, b, t);
        Vector3 p1 = Vector3.Lerp(b, c, t);
        return Vector3.Lerp(p0, p1, t);
    }

    public static Vector3 EvaluateCubic(Vector3 a, Vector3 b, Vector3 c, Vector3 d, float t)
    {
        Vector3 p0 = EvaluateQuadratic(a, b, c, t);
        Vector3 p1 = EvaluateQuadratic(b, c, d, t);
        return Vector3.Lerp(p0, p1, t);
    }

    public static Vector3 Tangent(Vector3 a, Vector3 b, Vector3 c, Vector3 d, float t)
    {
        return 3f * (1 - t) * (1 - t) * (b - a) + 6f * t * (1 - t) * (c - b) + 3f * t * t * (d - c);
    }

    private static Vector2 EvaluateQuadratic(Vector2 a, Vector2 b, Vector2 c, float t)
    {
        Vector3 p0 = Vector3.Lerp(a, b, t);
        Vector3 p1 = Vector3.Lerp(b, c, t);
        return Vector3.Lerp(p0, p1, t);
    }

    public static Vector3 EvaluateCubic(Vector2 a, Vector2 b, Vector2 c, Vector2 d, float t)
    {
        Vector3 p0 = EvaluateQuadratic(a, b, c, t);
        Vector3 p1 = EvaluateQuadratic(b, c, d, t);
        return Vector3.Lerp(p0, p1, t);
    }

    public static Vector3 Tangent(Vector2 a, Vector2 b, Vector2 c, Vector2 d, float t)
    {
        return 3f * (1 - t) * (1 - t) * (b - a) + 6f * t * (1 - t) * (c - b) + 3f * t * t * (d - c);
    }
}
