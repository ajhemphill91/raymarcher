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

// l < sin(a) < l + u
float sinp(float a, float l, float u)
{
    return l+((u/2.0)*(sin(a)+u));
}

//Change between integer values every t/a seconds, a different possible vals
float domains(float t, float a)
{
    const float PI = 3.1415926;
    return floor(sinp(2.0*PI/t*iGlobalTime, 0.0, a));
}

float rubiksCube(vec3 p)
{
    // Slowly spin the whole scene
    float PI = 3.1415926;
    float f = 2.0*PI*iGlobalTime*0.3;
    float phase = PI*0.3;

    //Active rotation should be switched here before the geom loop
    mat3 appliedRot, invAppliedRot;

    appliedRot = IdMatrix();
    invAppliedRot = IdMatrix();

    p = Rz(iGlobalTime / 2.0) * Rx(iGlobalTime / 2.0) * p;
    //p = Ry(-PI/2.0)*p;
    float r = 1.0;
    float space = 0.5;//(sin(f)+1.0)/2.0 + 0.05;
    float d = 2.0*r + space;
    float boxes;
    vec3 o;
    vec3 dim = vec3(r, r, r);

    //Timed rotations
    float dom = domains(4.0, 0.0);
    int which = 0;
    if(dom == 0.0)
    {
        which = 0;
    }
    else if(dom == 1.0)
    {
        which = 1;
    }
    else if(dom == 2.0)
    {
        which = 2;
    }


    for(int j = 0; j < 3; j++)
    {
      if(which == 0)
      {
        if(j == 0)
        {
          appliedRot = Ry(f);
          invAppliedRot = Ry(-f);
          o.y = -d;
        }
        else
        {
          appliedRot = IdMatrix();
          invAppliedRot = IdMatrix();
        }
      }
      else
      {
        if(j == 0)
        {
          //groupRotY = yCCW;
          o.y = -d;
        }
      }

      if(j == 1)
      {
        //groupRotY = IdMatrix();
        o.y = 0.0;
      }
      if(j == 2)
      {
        //groupRotY = yCW;
        //groupRotY = IdMatrix();
        o.y = d;
      }

      for(int k = 0; k < 3; k++)
      {
        if(which == 1)
        {
          if(k == 0)
          {
              appliedRot = Rz(-f);
              invAppliedRot = Rz(f);
              o.z = d;
          }
          else
          {
              appliedRot = IdMatrix();
              invAppliedRot = IdMatrix();
          }
        }
        else
        {
          if(k == 0)
          {
            //groupRotZ = zCW;
            //groupRotZ = IdMatrix();
            o.z = d;
          }
        }

        if(k == 1)
        {
          //groupRotZ = IdMatrix();
          o.z = 0.0;
        }
        if(k == 2)
        {
          //groupRotZ = zCCW;
          o.z = -d;
        }
        for(int i = 0; i < 3; i++)
        {
          //The specific rotation and it's inverse case should be assigned here for (i,j,k)
          if(which == 2)
          {
            if(i == 0)
            {
              //groupRotX = xCCW;
              o.x = -d;
              appliedRot = Rx(f);
              invAppliedRot = Rx(-f);
            }
            else
            {
              appliedRot = IdMatrix();
              invAppliedRot = IdMatrix();
            }
          }
          else
          {
            if(i == 0)
            {
              o.x = -d;
            }
          }

          if(i == 1)
          {
            //groupRotX = IdMatrix();
            o.x = 0.0;
          }
          if(i == 2)
          {
            //groupRotX = xCW;
            //groupRotX = IdMatrix();
            o.x = d;
          }

          //Render a cube at (i,j,k) - R*o gives the column/row rotation, and the invR*(p+R*o) keeps the cubes aligned with eachother during rotation
          if(i == 0 && j == 0 && k == 0)
          {
            boxes = roundBox(invAppliedRot*(p + appliedRot*o), dim, 0.1);
          }
          else
          {
            boxes = unionSDF(boxes, roundBox(invAppliedRot*(p + appliedRot*o), dim, 0.1));
          }

        }
      }
    }

    return boxes;
}

float sceneSDF(vec3 p)
{
    return rubiksCube(p);
}

/**
    Raymarching engine
*/

float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end)
{
    float depth = start;
    const int MAX_MARCHING_STEPS = 100;
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
    vec3 center = vec3(-5.0, -5.0, 0.0);
    vec3 up = vec3(0.0, 1.0, 0.0);
    mat3 viewToWorld = viewMat(eye, center, up);
    vec3 worldDir = viewToWorld * viewDir;

    //Cast the ray into the scene (entrypoint to scene())
    const float MIN_DIST = 0.0;
    const float MAX_DIST = 50.0;
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
