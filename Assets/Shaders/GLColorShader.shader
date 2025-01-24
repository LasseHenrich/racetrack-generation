Shader "Custom/GLColorShader"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        Pass
        {
            ZWrite Off
            Lighting Off
            Cull Off
            Blend SrcAlpha OneMinusSrcAlpha

            SetTexture[_MainTex] { combine primary }
        }
    }
}
