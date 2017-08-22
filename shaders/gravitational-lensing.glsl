/**
    3D rotation matrices
*/

mat3 Rx(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(1, 0, 0),
        vec3(0, c, -s),
        vec3(0, s, c)
    );
}

mat3 Ry(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

mat3 Rz(float theta)
{
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, -s, 0),
        vec3(s, c, 0),
        vec3(0, 0, 1)
    );
}

mat3 IdMatrix()
{
    return mat3(
            vec3(1, 0, 0),
            vec3(0, 1, 0),
            vec3(0, 0, 1)
    );
}

/**
    Signed distance surface combination functions
*/

float intersectSDF(float a, float b)
{
    return max(a, b);
}

float unionSDF(float a, float b)
{
    return min(a, b);
}

float diffSDF(float a, float b)
{
    return max(a, -b);
}

/**
    Signed distance geometry functions
*/

float box(vec3 p, vec3 size)
{
    vec3 d = abs(p) - (size / 2.0);
    float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
    float outsideDistance = length(max(d, 0.0));

    return insideDistance + outsideDistance;
}

float roundBox(vec3 p, vec3 halfFaceWidth, float r)
{
    return length(max(abs(p)-halfFaceWidth,0.0))-r;
}

float cylinder(vec3 p, float h, float r)
{
    float inOutRadius = length(p.xy) - r;
    float inOutHeight = abs(p.z) - h/2.0;
    float insideDistance = min(max(inOutRadius, inOutHeight), 0.0);
    float outsideDistance = length(max(vec2(inOutRadius, inOutHeight), 0.0));

    return insideDistance + outsideDistance;
}

float sdPlane( vec3 p, vec4 n )
{
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}

float sdSphere( vec3 p, float s )
{
  vec3 a = sin(p)/mod(p.z, s);
  vec3 b = cos(p)/mod(p.x, s);
  vec3 c = sin(p)/mod(p.y, s);
  return length(p)-s-length(dot(cross(a,b),c));
}

float cylinder2(vec3 p, float h, float r)
{
    float PI = 3.14159;
    p = mod(cos(p)*2.0-1.0, vec3(32.0*PI));
    p.z = 0.5*(sin(iGlobalTime)*p.z*2.0 - 1.0);

    float ior = length(p.xy) - r;
    float ioh = abs(p.z) - h/2.0;
    float id = min(max(ior, ioh), 10.0);
    float od = length(max(vec2(ior, ioh), 0.0));
    return od + id;
}


float meltingCylinder(vec3 p)
{
    float PI = 3.14159;
    p = Rz(PI/2.)*Ry(PI/2.)*Rx(sin(iGlobalTime))*p;
    float h = 200.0;
    float r = 0.25;
    float cyl = cylinder2(p, h, r);
    float s = 2.0;
    float sphere = sdSphere(p,s);
    float res = unionSDF(cyl, sphere);

    return res;
}

float sceneSDF(vec3 p)
{
    return meltingCylinder(p);
}

/**
    Raymarching engine
*/

float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end)
{
    float depth = start;
    const int MAX_MARCHING_STEPS = 255;
    float EPSILON = 0.0001;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++)
    {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON)
        {
			      return depth;
        }
        depth += dist;
        if (depth >= end)
        {
            return end;
        }
    }
    return end;
}

vec3 rayDirection(float fov, vec2 size, vec2 uv)
{
    //uv = uv*2.0 - 1.0;
    vec2 xy = uv - size / 2.0;
    float z = size.y / tan(radians(fov) / 2.0);
    return normalize(vec3(xy, -z));
}

/**
    Lighting functions
*/

vec3 estimateNormal(vec3 p)
{
    float EPSILON = 0.0001;
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

vec3 phongContribForLight(vec3 k_d, vec3 k_s, float shininess, vec3 p,
                          vec3 eye, vec3 lightPos, vec3 lightIntensity)
{
    vec3 N = estimateNormal(p);
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));

    float dotLN = dot(L, N);
    float dotRV = dot(R, V);

    if (dotLN < 0.0) return vec3(0.0, 0.0, 0.0);
    if (dotRV < 0.0) return lightIntensity * (k_d * dotLN);
    return lightIntensity * ( (k_d * dotLN) + (k_s * pow(dotRV, shininess)) );
}

vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float shininess, vec3 p, vec3 eye)
{
    vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;

    vec3 light1Pos = vec3(0.0, 15.0, 0.0);
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);

    color += phongContribForLight(k_d, k_s, shininess, p, eye,
                                  light1Pos,
                                  light1Intensity);

    return color;
}


mat3 viewMat(vec3 eye, vec3 center, vec3 up)
{
    vec3 f = normalize(center - eye);
    vec3 s = normalize(cross(f, up));
    vec3 u = cross(s, f);
    return mat3(s, u, -f);
}

void main()
{
    //Compute a ray for this fragcoord in camera coords, xform to world coords
    vec3 viewDir = rayDirection(45.0, iResolution.xy, gl_FragCoord.st);
    vec3 eye = vec3(0.0, 0.0, 30.0);
    vec3 center = vec3(0.0, 0.0, 0.0);
    vec3 up = vec3(0.0, 1.0, 0.0);
    mat3 viewToWorld = viewMat(eye, center, up);
    vec3 worldDir = viewToWorld * viewDir;

    //Cast the ray into the scene (entrypoint to scene())
    const float MIN_DIST = 0.0;
    const float MAX_DIST = 100.0;
    float dist = shortestDistanceToSurface(eye, worldDir, MIN_DIST, MAX_DIST);

    //No intersections
    const float EPSILON = 0.0001;
    if (dist > MAX_DIST - EPSILON)
    {
          gl_FragColor = vec4(0.0, 0.0, 0.0, 0.0);
		      return;
    }

    //Closest point p of the signed distance function
    vec3 p = eye + dist * worldDir;

    // Use the surface normal as the ambient color of the material
    vec3 K_a = (estimateNormal(p) + vec3(1.0)) / 2.0;
    vec3 K_d = K_a;
    vec3 K_s = vec3(1.0, 1.0, 1.0);
    float shininess = 10.0;

    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, eye);

    gl_FragColor = vec4(color, 1.0);

}
