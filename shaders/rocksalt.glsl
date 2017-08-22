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

float smoothUnionSDF(float a, float b, float k)
{
    float res = exp(-k*a)+exp(-k*b);
    return -log(res)/k;
}

float polynomialUnionSDF(float a, float b, float k)
{
    float h = clamp(0.5+0.5*(b-a)/k, 0.0, 1.0);
    return mix(b, a, h) - k*h*(1.0-h);
}

float powerUnionSDF(float a, float b, float k)
{
    a = pow(a, k);
    b = pow(b, k);
    return pow( (a*b)/(a+b), 1.0/k);
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

vec3 opRep(vec3 p, vec3 c)
{
    return mod(p,c)-0.5*c;
}

float sceneSDF(vec3 p)
{
    float scene;
    const float PI = 3.1415926;
    p = Rx(PI/2.0)*p;
    vec3 c = vec3(9.0, 9.0, 9.0);
    p = opRep(p, c);
    float boxes = roundBox(p, vec3(0.5, 0.5, 0.5), 0.2);


    float pillarHeight = 10.0;
    float pillarRadius = 0.1;
    float pillars = cylinder(p, pillarHeight, pillarRadius);
    p = Rx(-PI/2.0)*p;
    pillars = polynomialUnionSDF(pillars, cylinder(p, pillarHeight, pillarRadius), 1.);
    p = Ry(PI/2.0)*p;
    pillars = smoothUnionSDF(pillars, cylinder(p, pillarHeight/2., pillarRadius), sin(iGlobalTime)*0.5+2.0);
    vec3 oldP = p;
    p = cross(opRep(p, c), opRep(c, p));
    pillars = polynomialUnionSDF(pillars, cylinder(p, pillarHeight, pillarRadius), sin(4.0*iGlobalTime)*0.1+1.);
    pillars = smoothUnionSDF(pillars, boxes, abs(sin(iGlobalTime))+1.0);
    return pillars;
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
    uv = uv*2.0 -1.0;
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

    vec3 light1Pos = vec3(15.0, 15.0, 20.0);
    vec3 light1Intensity = vec3(0.1, 0.1, 0.1);

    color += phongContribForLight(k_d, k_s, shininess, p, eye,
                                  eye+light1Pos,
                                  light1Intensity);

    vec3 light2Pos = vec3(-15.0, -15.0, 50.0);
    vec3 light2Intensity = vec3(0.1, 0.1, 0.1);

    color += phongContribForLight(k_d, k_s, shininess, p, eye,
                                  eye+light2Pos,
                                  light2Intensity);

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
    vec3 eye = vec3(iGlobalTime, 0.0, iGlobalTime+2.);
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
    vec3 N = estimateNormal(p);
    vec3 K_a = (estimateNormal(p) + vec3(1.0)) / 2.0;
    vec3 K_d = K_a;
    float NX = dot(N.xz, gl_FragCoord.st);
    float NY = dot(N.yz, gl_FragCoord.st);
    //vec3 K_a = vec3(0.0, 0.0, 0.0);
    //vec3 K_d = vec3(0.6, 1.0-NY, 1.0-NX);
    vec3 K_s = vec3(0.2, 0.2, 0.2);
    float shininess = 1.0;

    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, eye);

    gl_FragColor = vec4(color, 1.0);

}
