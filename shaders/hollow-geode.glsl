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
    Euclidian norm with higher exponents
*/
float length8(vec2 p)
{
    p = p*p; p = p*p; p = p*p;
	  return pow(p.x + p.y, 1.0/8.0);
}

float length8(vec3 p)
{
    p = p*p; p = p*p; p = p*p;
    return pow(p.x + p.y + p.z, 1.0/8.0);
}

float length4(vec2 p)
{
    p = p*p; p = p*p;
    return pow(p.x + p.y, 1.0/4.0);
}

float length4(vec3 p)
{
    p = p*p; p = p*p;
    return pow(p.x + p.y + p.z, 1.0/8.0);
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

float diffSDF(float a, float b)
{
    return max(a, -b);
}

/*
float opScale(vec3 p, float s)
{
    return primitive(p/s)*s;
}
*/

vec3 opRep(vec3 p, vec3 c)
{
    return mod(p, c)-0.5*c;
}

float sinDisplace(vec3 p)
{
    return (sin(iGlobalTime*p.x)+sin(iGlobalTime*p.y)+sin(iGlobalTime*p.z));
}

vec3 opTwist(vec3 p, float amt)
{
    float c = cos(amt*p.y);
    float s = sin(amt*p.y);
    mat2  m = mat2(c,-s,s,c);
    return vec3(m*p.xz, p.y);
}

vec3 opBend(vec3 p, float amt)
{
    float c = cos(amt*p.y);
    float s = sin(amt*p.y);
    mat2  m = mat2(c,-s,s,c);
    return vec3(m*p.xy, p.z);
}

/**
    Signed distance geometry functions
*/

float sphere(vec3 p, float s)
{
    return length(p)-s;
}

float box(vec3 p, vec3 b)
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float roundBox(vec3 p, vec3 halfFaceWidth, float cornerRadius)
{
    return length(max(abs(p)-halfFaceWidth,0.0))-cornerRadius;
}

float infCylinder(vec3 p, vec3 c)
{
    return length(p.xz-c.xy)-c.z;
}

float cylinder(vec3 p, vec2 rh)
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - rh;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float cone(vec3 p, vec2 c)
{
    // c must be normalized
    float q = length(p.xy);
    return dot(c,vec2(q,p.z));
}

float plane(vec3 p, vec4 n)
{
    // n must be normalized
    return dot(p,n.xyz) + n.w;
}

float ellipsoid(vec3 p, vec3 r)
{
    return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z);
}

float triPrism(vec3 p, vec2 h)
{
    vec3 q = abs(p);
    return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
}

float hexPrism(vec3 p, vec2 h)
{
    vec3 q = abs(p);
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
}

//Capsules, a and b are the endpoints
float capsule(vec3 p, vec3 a, vec3 b, float r)
{
    vec3 pa = p - a;
    vec3 ba = b - a;
    float h = clamp(dot(pa,ba)/dot(ba,ba), 0.0, 1.0);
    return length(pa - ba*h) - r;
}

//Torus: t controls width, depth
float torus(vec3 p, vec2 t)
{
    vec2 q = vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}

float torus82(vec3 p, vec2 t)
{
    vec2 q = vec2(length8(p.xz)-t.x, p.y);
    return length(q)-t.y;
}

float torus88(vec3 p, vec2 t)
{
    vec2 q = vec2(length8(p.xz)-t.x, p.y);
    return length8(q)-t.y;
}

/**
    Useful sinusoids
*/

// l < sin(a) < l + u
float sinp(float a, float l, float u)
{
    return l+((u/2.0)*(sin(a)+u));
}

float sint(float l, float u)
{
    return sinp(iGlobalTime, l, u);
}

// l < cos(a) < l + u
float cosp(float a, float l, float u)
{
    return l+((u/2.0)*(cos(a)+u));
}

float cost(float l, float u)
{
    return cosp(iGlobalTime, l, u);
}

//Change between integer values every t/a seconds, a different possible vals
float domains(float t, float a)
{
    const float PI = 3.1415926;
    return floor(sinp(2.0*PI/t*iGlobalTime, 0.0, a));
}

vec3 simplexNoise()
{
    vec3 z = vec3(0.0);
    vec2 uv = gl_FragCoord.st/iResolution.xy;
    float r = 115.47*uv.x;
    vec2 s = fract(vec2(r, uv.y + r*0.5));
    if(s.x > s.y)
    {
      z.xy = 1.0-vec2(s.x, s.y-s.x);
      z.z = s.y;
    }
    else
    {
      z.yz = 1.0-vec2(s.x-s.y, s.y);
      z.x = s.x;
    }

    return fract(z);
}

/**
    The scene
*/

float sceneSDF(vec3 p)
{
    float scene = 0.0;
    const float PI = 3.1415926;

    p.x -= 2.5;
    p.y -= 2.5;
    p = Rx(iGlobalTime/2.)*Ry(iGlobalTime/2.)*p;
    p = Rx(PI/3.0)*p;

    float wall = roundBox(p, vec3(5., 5., 2.0), 0.5);
    vec3 p0 = opTwist(p, 0.1);
    float spheres = sphere(vec3(p0.x+5., p0.y+5., p0.z), 2.5);
    spheres = unionSDF(spheres, sphere(vec3(p0.x-5., p0.y+5., p0.z), 2.5));
    spheres = unionSDF(spheres, sphere(vec3(p0.x+5., p0.y-5., p0.z), 2.5));
    spheres = unionSDF(spheres, sphere(vec3(p0.x-5., p0.y-5., p0.z), 2.5));
    wall = intersectSDF(wall, spheres);

    float nuc = hexPrism(vec3(p.x, p.y, p.z+2.1), vec2(1.25, 0.1));
    nuc = unionSDF(nuc, hexPrism(vec3(p.x, p.y, p.z-2.1), vec2(1.25, 0.1)));
    p = Ry(PI/2.)*p;
    nuc = unionSDF(nuc, hexPrism(vec3(p.x, p.y, p.z+5.25), vec2(0.5, 0.1)));
    nuc = unionSDF(nuc, hexPrism(vec3(p.x, p.y, p.z-5.25), vec2(0.5, 0.1)));
    p = Ry(PI/-2.)*p;
    wall = polynomialUnionSDF(wall, nuc, 7.);
    float cx = dot(p.y, p.z);
    float cy = dot(p.x, p.z);
    float cz = dot(p.x, p.y);
    vec3 c = vec3(sint(cx, cz*cx), sint(cy, cz/cy), sint(cz, 0.25*cz));
    p = opRep(p, c);
    float s = PI/3.;
    float doors = ellipsoid(p, vec3(s, s, 2.*s));

    scene = diffSDF(wall, doors);

    return scene;
}

/**
    Raymarching engine
*/

float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end)
{
    float depth = start;
    const int MAX_MARCHING_STEPS = 255;
    const float EPSILON = 0.0001;
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
    const float EPSILON = 0.0001;
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

    vec3 light1Pos = vec3(eye.x, eye.y, eye.z);
    vec3 light1Intensity = vec3(0.6, 0.6, 0.6);

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

vec3 unpackColor(float f)
{
    vec3 color;
    color.b = floor(f / 256.0 / 256.0);
    color.g = floor((f - color.b * 256.0 * 256.0) / 256.0);
    color.r = floor(f - color.b * 256.0 * 256.0 - color.g * 256.0);
    // now we have a vec3 with the 3 components in range [0..255]. Let's normalize it!
    return color / 255.0;
}

void main()
{
    //Compute a ray for this fragcoord in camera coords, xform to world coords
    vec2 uv = gl_FragCoord.st/iResolution.xy;
    vec3 viewDir = rayDirection(45.0, iResolution.xy, gl_FragCoord.st);
    vec3 eye = vec3(0.0, 0.0, 25.);
    vec3 center = vec3(-1.75, -2.0, 0.0);
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
          vec3 bg = fract(simplexNoise()*iGlobalTime)/6.;
          gl_FragColor = vec4(0.0, uv.y-bg.g, uv.x-bg.b, 0.0);
		      return;
    }

    //Closest point p of the signed distance function
    vec3 p = eye + dist * worldDir;
    // Use the surface normal as the ambient color of the material
    vec3 n = simplexNoise()*(cos(uv.x) + sin(uv.y));
    vec3 K_a = (estimateNormal(p) + vec3(1.0)) / 6.;
    vec3 K_d = fract(cross(p, cross(n, -n)))/2.;
    vec3 K_s = vec3(0.75, 0.75, 0.75);
    float shininess = 500.0;

    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, eye);
    //vec3 color = unpackColor(dist*256.);

    gl_FragColor = vec4(color, 1.0);

}
