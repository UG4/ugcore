/*
 * neurite_projector.cpp
 *
 *  Created on: 19.12.2016
 *      Author: mbreit
 */

#include "neurite_projector.h"
#include "common/error.h"
#include "lib_grid/global_attachments.h"

#include <boost/lexical_cast.hpp>

namespace ug {


NeuriteProjector::NeuriteProjector()
: m_quadOrder(80)
{
    // If branching points extend 5r in each direction and we integrate over twice that size
    // we add around quadOrder/2 Gaussians with sigma=r over a total length of 20r.
    // So we should use about quadOrder=80 to ensure a smooth surface
    // that does not look like a pearl necklace.
    prepare_quadrature();
}


NeuriteProjector::NeuriteProjector(SPIGeometry3d geometry)
: RefinementProjector(geometry),
  m_quadOrder(80)
{
    attach_surf_params();
    prepare_quadrature();
}


NeuriteProjector::~NeuriteProjector()
{}


void NeuriteProjector::set_geometry(SPIGeometry3d geometry)
{
    // call base class method
    RefinementProjector::set_geometry(geometry);

    attach_surf_params();
}



number NeuriteProjector::new_vertex(Vertex* vrt, Edge* parent)
{
    return push_onto_surface(vrt, parent);
}


number NeuriteProjector::new_vertex(Vertex* vrt, Face* parent)
{
    return push_onto_surface(vrt, parent);
}


void NeuriteProjector::direction_at_grid_object(vector3& dirOut, GridObject* o) const
{
    // treat vertex separately as it is no vertex group
    if (o->base_object_id() == VERTEX)
    {
        Vertex* v = dynamic_cast<Vertex*>(o);
        UG_COND_THROW(!v, "Non-vertex with VERTEX base object id.")

        const SurfaceParams& sp = m_aaSurfParams[v];
        uint32 nid = sp.neuriteID;
        float t = sp.axial;

        std::vector<Section>::const_iterator secIt = get_section_it(nid, t);

        number te = secIt->endParam;
        const number* s = &secIt->splineParamsX[0];
        number& v0 = dirOut[0];
        v0 = -3.0*s[0]*(te-t) - 2.0*s[1];
        v0 = v0*(te-t) - s[2];

        s = &secIt->splineParamsY[0];
        number& v1 = dirOut[1];
        v1 = -3.0*s[0]*(te-t) - 2.0*s[1];
        v1 = v1*(te-t) - s[2];

        s = &secIt->splineParamsZ[0];
        number& v2 = dirOut[2];
        v2 = -3.0*s[0]*(te-t) - 2.0*s[1];
        v2 = v2*(te-t) - s[2];

        return;
    }

    IVertexGroup* vrtGrp = dynamic_cast<IVertexGroup*>(o);
    UG_COND_THROW(!vrtGrp, "Non-vertex element which is not a vertex group.");

    size_t nid;
    float t, dummy;
    average_params(nid, t, dummy, vrtGrp);

    std::vector<Section>::const_iterator secIt = get_section_it((uint32)nid, t);

    number te = secIt->endParam;
    const number* s = &secIt->splineParamsX[0];
    number& v0 = dirOut[0];
    v0 = -3.0*s[0]*(te-t) - 2.0*s[1];
    v0 = v0*(te-t) - s[2];

    s = &secIt->splineParamsY[0];
    number& v1 = dirOut[1];
    v1 = -3.0*s[0]*(te-t) - 2.0*s[1];
    v1 = v1*(te-t) - s[2];

    s = &secIt->splineParamsZ[0];
    number& v2 = dirOut[2];
    v2 = -3.0*s[0]*(te-t) - 2.0*s[1];
    v2 = v2*(te-t) - s[2];
}



void NeuriteProjector::attach_surf_params()
{
    Grid& grid = this->geometry()->grid();

    UG_COND_THROW(!GlobalAttachments::is_declared("npSurfParams"),
        "GlobalAttachment 'npSurfParams' not declared.");
    m_aSurfParams = GlobalAttachments::attachment<Attachment<SurfaceParams> >("npSurfParams");

    // make sure surfaceParams attachment is attached
    UG_COND_THROW(!grid.has_vertex_attachment(m_aSurfParams),
        "No surface parameter attachment for neurite projector attached to grid.");

    m_aaSurfParams.access(grid, m_aSurfParams);

    // handle attachment values also on higher grid levels (if required)
    MultiGrid* mg = dynamic_cast<MultiGrid*>(&grid);
    if (mg)
    {
        SmartPtr<MultiGrid> spMG(mg);

        // never destroy the grid from here - we did not create it
        ++(*spMG.refcount_ptr());

        m_cah.set_attachment(m_aSurfParams);
        m_cah.set_grid(spMG);
    }
}



void NeuriteProjector::add_neurite(const Neurite& n)
{
    m_vNeurites.push_back(n);
}


const NeuriteProjector::Neurite& NeuriteProjector::neurite(size_t nid) const
{
    UG_COND_THROW(nid >= m_vNeurites.size(),
        "Requested neurite index " << nid << " that is not present in the neurite.");
    return m_vNeurites[nid];
}


Grid::VertexAttachmentAccessor<Attachment<NeuriteProjector::SurfaceParams> >&
NeuriteProjector::surface_params_accessor()
{
    return m_aaSurfParams;
}


const std::vector<std::pair<number, number> >& NeuriteProjector::quadrature_points() const
{
    return m_qPoints;
};



std::vector<NeuriteProjector::Section>::const_iterator
NeuriteProjector::get_section_it(uint32 nid, float t) const
{
    Section cmpSec(t);
    const std::vector<Section>& vSections = m_vNeurites[nid].vSec;
    std::vector<Section>::const_iterator itSec =
        std::lower_bound(vSections.begin(), vSections.end(), cmpSec, CompareSections());

    UG_COND_THROW(itSec == vSections.end(),
        "Could not find section for parameter t = " << t << " in neurite " << nid << ".");

    return itSec;
}



static bool cmpQPairs(const std::pair<number, number>& a, const std::pair<number, number>& b)
{return a.first < b.first;}

void NeuriteProjector::prepare_quadrature()
{
    GaussLegendre gl(m_quadOrder);
    size_t nPts = gl.size();
    m_qPoints.resize(nPts);
    for (size_t i = 0; i < nPts; ++i)
    {
        m_qPoints[i].first = gl.point(i)[0];
        m_qPoints[i].second = gl.weight(i);
    }

    std::sort(m_qPoints.begin(), m_qPoints.end(), cmpQPairs);
}


void NeuriteProjector::average_params(size_t& neuriteID, float& t, float& angle, const IVertexGroup* parent) const
{
    // branching vertices need to belong to parent branch
    // AND parent branch neurite index must be smaller than child's index
    neuriteID = 0;
    t = 0.0;
    vector2 v(0.0);
    size_t nVrt = parent->num_vertices();
    for (size_t i = 0; i < nVrt; ++i)
    {
        const SurfaceParams& surfParams = m_aaSurfParams[parent->vertex(i)];
        size_t nid = (size_t) surfParams.neuriteID;
        neuriteID = std::max(neuriteID, nid);
    }
    for (size_t i = 0; i < nVrt; ++i)
    {
        const SurfaceParams& surfParams = m_aaSurfParams[parent->vertex(i)];
        size_t nid = (size_t) surfParams.neuriteID;
        if (nid == neuriteID)
        {
            t += surfParams.axial;
            float phi = surfParams.angular;
            VecAdd(v, v, vector2(cos(phi), sin(phi)));
        }
    }
    t /= nVrt;

    if (fabs(v[0]) < 1e-8)
    {
      if (fabs(v[1]) < 1e-8)
          UG_THROW("Angle for new vertex cannot be determined.")
      else
          angle = v[1] < 0 ? 1.5*PI : 0.5*PI;
    }
    else
       angle = v[0] < 0 ? PI - atan(-v[1]/v[0]) : atan(v[1]/v[0]);

    if (angle < 0) angle += 2.0*PI;
}


void NeuriteProjector::average_pos_from_parent(vector3& posOut, const IVertexGroup* parent) const
{
    posOut = 0.0;
    size_t nVrt = parent->num_vertices();
    for (size_t i = 0; i < nVrt; ++i)
        posOut += this->pos(parent->vertex(i));
    posOut /= (number) nVrt;
}



static void compute_position_and_velocity_in_section
(
    vector3& posAxOut,
    vector3& velOut,
    number& radiusOut,
    std::vector<NeuriteProjector::Section>::const_iterator secIt,
    float t
)
{
    number te = secIt->endParam;
    const number* s = &secIt->splineParamsX[0];
    number& p0 = posAxOut[0];
    number& v0 = velOut[0];
    p0 = s[0]*(te-t) + s[1];  v0 = -3.0*s[0]*(te-t) - 2.0*s[1];
    p0 = p0*(te-t) + s[2];    v0 = v0*(te-t) - s[2];
    p0 = p0*(te-t) + s[3];

    s = &secIt->splineParamsY[0];
    number& p1 = posAxOut[1];
    number& v1 = velOut[1];
    p1 = s[0]*(te-t) + s[1];  v1 = -3.0*s[0]*(te-t) - 2.0*s[1];
    p1 = p1*(te-t) + s[2];    v1 = v1*(te-t) - s[2];
    p1 = p1*(te-t) + s[3];

    s = &secIt->splineParamsZ[0];
    number& p2 = posAxOut[2];
    number& v2 = velOut[2];
    p2 = s[0]*(te-t) + s[1];  v2 = -3.0*s[0]*(te-t) - 2.0*s[1];
    p2 = p2*(te-t) + s[2];    v2 = v2*(te-t) - s[2];
    p2 = p2*(te-t) + s[3];

    s = &secIt->splineParamsR[0];
    radiusOut = s[0]*(te-t) + s[1];
    radiusOut = radiusOut*(te-t) + s[2];
    radiusOut = radiusOut*(te-t) + s[3];
}


static void bp_defect_and_gradient
(
    number& defectOut,
    vector3& gradientOut,
    const std::vector<NeuriteProjector::BPProjectionHelper>& bpList,
    const vector3& x,
    const NeuriteProjector* np
)
{
    // get integration rule points and weights
    const std::vector<std::pair<number, number> >& qPoints = np->quadrature_points();
    size_t nPts = qPoints.size();

    // evaluate integrand at points
    defectOut = 0.0;
    gradientOut = 0.0;
    size_t nParts = bpList.size();
    for (size_t i = 0; i < nParts; ++i)
    {
        const NeuriteProjector::BPProjectionHelper& helper = bpList[i];
        const number& start = helper.start;
        const number& end = helper.end;
        std::vector<NeuriteProjector::Section>::const_iterator secIt = helper.sec_start;

        // iterate IPs
        for (size_t j = 0; j < nPts; ++j)
        {
            // transform point coords to current context
            float t = (float) (start + qPoints[j].first*(end-start));
            while (secIt->endParam < t) ++secIt;    // this should be safe as last section must end with 1.0

            vector3 posAx, vel;
            number radius;
            compute_position_and_velocity_in_section(posAx, vel, radius, secIt, t);

            number radInv = 1.0 / radius;
            VecSubtract(posAx, posAx, x);
            number posNormSq = VecNormSquared(posAx);
            number velNorm = sqrt(VecNormSquared(vel));
            number integrandVal = radInv*exp(-posNormSq*radInv*radInv) * velNorm;
            number gradIntegrandVal = 2.0*radInv*radInv * integrandVal;
            VecScale(posAx, posAx, gradIntegrandVal);

            number w = qPoints[j].second;
            defectOut += integrandVal * w * (end-start);
            VecScaleAdd(gradientOut, 1.0, gradientOut, w * (end-start), posAx);

            // debugging
            //UG_LOGN("  Ival: " << integrandVal << ",  iExp: " << -posNormSq*radInv*radInv
            //    << ",  velNorm: " << velNorm << ",  weight: " << w
            //    << ",  t: " << t << ",  intvl: " << end-start);
        }
    }

    defectOut -= sqrt(PI)*exp(-1.0);
}


static inline void compute_ONB
(
    vector3& firstOut,
    vector3& secondOut,
    vector3& thirdOut,
    const vector3& firstIn,
    const vector3& refDir
)
{
    VecNormalize(firstOut, firstIn);
    number fac = VecProd(refDir, firstOut);
    VecScaleAdd(secondOut, 1.0, refDir, -fac, firstOut);
    VecNormalize(secondOut, secondOut);

    VecCross(thirdOut, firstOut, secondOut);
}


// computes the polar angle phi of a vector
// given a rotational axis direction and two additional
// directions, together being a orthoNORMAL basis
static inline number compute_angle
(
    const vector3& axis,
    const vector3& secondAxis,
    const vector3& thirdAxis,
    vector3& posFromAxis
)
{
    // eliminate component in axis direction
    VecScaleAdd(posFromAxis, 1.0, posFromAxis, -VecProd(posFromAxis, axis), axis);

    // compute coord vector (components 2 and 3) w.r.t. ONB
    vector2 relCoord;
    relCoord[0] = VecProd(posFromAxis, secondAxis);
    VecScaleAdd(posFromAxis, 1.0, posFromAxis, -relCoord[0], secondAxis);
    relCoord[1] = VecProd(posFromAxis, thirdAxis);
    VecNormalize(relCoord, relCoord);

    // get angle from coord vector
    number angle;
    if (fabs(relCoord[0]) < 1e-8)
        angle = relCoord[1] < 0 ? 1.5*PI : 0.5*PI;
    else
        angle = relCoord[0] < 0 ? PI - atan(-relCoord[1]/relCoord[0]) : atan(relCoord[1]/relCoord[0]);

    // angle should be in [0,2*PI)
    if (angle < 0) angle += 2.0*PI;

    return angle;
}


static void newton_for_bp_projection
(
    vector3& posOut,
    const std::vector<NeuriteProjector::BPProjectionHelper>& vProjHelp,
    const NeuriteProjector* np
)
{
    // calculate start defect and gradient
    number def;
    number def_init;
    vector3 grad;
    bp_defect_and_gradient(def, grad, vProjHelp, posOut, np);
    def_init = fabs(def);

    // perform some kind of Newton search for the correct position
    size_t maxIt = 100;
    number minDef = 1e-8*sqrt(PI)*exp(-1.0);
    size_t iter = 0;
    while (fabs(def) > minDef && fabs(def) > 1e-8 * def_init && ++iter <= maxIt)
    {
//UG_LOGN(iter << "  " << def << "  " << grad << "  " << pos);
        vector3 posTest;
        vector3 gradTest;
        number defTest;
        VecScaleAdd(posTest, 1.0, posOut, -def / VecNormSquared(grad), grad);
        bp_defect_and_gradient(defTest, gradTest, vProjHelp, posTest, np);

        // line search
        number lambda = 0.5;
        number bestDef = defTest;
        vector3 bestGrad = gradTest;
        vector3 bestPos = posTest;
        while (fabs(defTest) >= fabs(def) && lambda > 0.001)
        {
//UG_LOGN("    line search with lambda = " << lambda);
            VecScaleAdd(posTest, 1.0, posOut, -lambda*def / VecNormSquared(grad), grad);
            bp_defect_and_gradient(defTest, gradTest, vProjHelp, posTest, np);
            if (fabs(defTest) < fabs(bestDef))
            {
                bestDef = defTest;
                bestGrad = gradTest;
                bestPos = posTest;
            }
            lambda *= 0.5;
        }
        def = bestDef;
        grad = gradTest;
        posOut = bestPos;
    }
    UG_COND_THROW(def != def,
        "Newton iteration did not converge for branching point projection. Defect is NaN.")
    UG_COND_THROW(fabs(def) > minDef && fabs(def) > 1e-8 * def_init,
        "Newton iteration did not converge for branching point projection.")
}



static void pos_on_surface_bp
(
    vector3& posOut,
    const NeuriteProjector::Neurite& neurite,
    float& t,
    float& angle,
    std::vector<NeuriteProjector::BranchingRegion>::const_iterator& it,
    const IVertexGroup* parent,
    const NeuriteProjector* np
)
{
    const std::vector<NeuriteProjector::Section>& vSections = neurite.vSec;

    // preparations for Newton method:
    // save integration start and end positions of all BRs of this BP,
    // also save a section iterator to start from
    SmartPtr<NeuriteProjector::BranchingPoint> bp = it->bp;
    size_t nParts = bp->vRegions.size();
    std::vector<NeuriteProjector::BPProjectionHelper> vProjHelp(nParts);
    const std::vector<NeuriteProjector::Section>* secs;
    for (size_t i = 0; i < nParts; ++i)
    {
        number& start = vProjHelp[i].start;
        number& end = vProjHelp[i].end;
        std::vector<NeuriteProjector::Section>::const_iterator& sec_start = vProjHelp[i].sec_start;

        size_t nid = bp->vNid[i];
        const NeuriteProjector::BranchingRegion* br = bp->vRegions[i];

        UG_COND_THROW(!br, "Requested branching region not accessible.");

        start = br->tstart;
        end = br->tend;
        secs = &np->neurite(nid).vSec;

        number dt = end - start;
        start = std::max(start - 0.5*dt, 0.0);
        end = std::min(end + 0.5*dt, 1.0);

        if (start == 0.0)
            sec_start = secs->begin();
        else
        {
            NeuriteProjector::Section cmpSec(start);
            sec_start = std::lower_bound(secs->begin(), secs->end(), cmpSec, NeuriteProjector::CompareSections());
            UG_COND_THROW(sec_start == secs->end(),
                "Could not find section for start parameter t = " << start << ".");
        }
    }

   // initial position: regular refinement
    np->average_pos_from_parent(posOut, parent);

    // perform Newton iteration
    try {newton_for_bp_projection(posOut, vProjHelp, np);}
    UG_CATCH_THROW("Newton iteration for neurite projection at branching point failed.");


    // update the surface param info to the new position
    // In order not to have to compute zeros of a 5th order polynomial,
    // we make a linearity approximation here:
    // s(t) = s(t0) + v(t0)*(t-t0);
    // v(t) = v(t0)
    // Then from v(t) * (s(t)-pos) = 0, we get:
    // t = t0 + v(t0) * (pos - s(t0)) / (v(t0)*v(t0))
    vector3 s, v;
    number dummy;
    std::vector<NeuriteProjector::Section>::const_iterator secIt = vSections.begin();
    while (t > secIt->endParam)
        ++secIt;
    compute_position_and_velocity_in_section(s, v, dummy, secIt, t);
    VecSubtract(s, posOut, s);
    t += VecProd(v,s)/VecProd(v,v);

    // and now angle
    vector3 projRefDir;
    vector3 thirdDir;
    compute_ONB(v, projRefDir, thirdDir, v, neurite.refDir);
    angle = compute_angle(v, projRefDir, thirdDir, s);

    // if we have ended up outside of the BP, then use the usual positioning
    if (t > it->tend)
    {
        vector3 posAx;
        number radius;
        compute_position_and_velocity_in_section(posAx, v, radius, secIt, t);
        VecScaleAdd(posOut, 1.0, posAx, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
    }
}


static void pos_on_surface_neurite
(
    vector3& posOut,
    const NeuriteProjector::Neurite& neurite,
    size_t neuriteID,
    float& t,
    float& angle
)
{
    const std::vector<NeuriteProjector::Section>& vSections = neurite.vSec;

    // find correct section
    NeuriteProjector::Section cmpSec(t);
    std::vector<NeuriteProjector::Section>::const_iterator it =
        std::lower_bound(vSections.begin(), vSections.end(), cmpSec, NeuriteProjector::CompareSections());

    UG_COND_THROW(it == vSections.end(),
        "Could not find section for parameter t = " << t << " in neurite " << neuriteID << ".");

    // calculate correct axial position
    // and velocity dir
    vector3 posAx, vel;
    number radius;
    compute_position_and_velocity_in_section(posAx, vel, radius, it, t);

    // calculate orthonormal basis vectors
    vector3 projRefDir, thirdDir;
    compute_ONB(vel, projRefDir, thirdDir, vel, neurite.refDir);

    // calculate new position
    VecScaleAdd(posOut, 1.0, posAx, radius*cos(angle), projRefDir, radius*sin(angle), thirdDir);
}


static void pos_on_surface_tip
(
    vector3& posOut,
    const NeuriteProjector::Neurite& neurite,
    const IVertexGroup* parent,
    const NeuriteProjector* np
)
{
    const std::vector<NeuriteProjector::Section>& vSections = neurite.vSec;

    // initial position: regular refinement
    np->average_pos_from_parent(posOut, parent);

    // project to half-sphere with given radius over tip center
    // TODO: One might think about something more sophisticated,
    // as the current approach ensures only continuity of the radius
    // at tips, but not differentiability.
    const NeuriteProjector::Section& sec = vSections[vSections.size()-1];
    vector3 tipCenter(sec.splineParamsX[3], sec.splineParamsY[3], sec.splineParamsZ[3]);
    number radius = sec.splineParamsR[3];

    vector3 radialVec;
    VecSubtract(radialVec, posOut, tipCenter);
    VecScale(radialVec, radialVec, radius/sqrt(VecProd(radialVec, radialVec)));
    VecAdd(posOut, tipCenter, radialVec);
}



number NeuriteProjector::push_onto_surface(Vertex* vrt, const IVertexGroup* parent)
{
    // average axial and angular params from parent;
    // also decide on neuriteID
    size_t neuriteID;
    float t;
    float angle;
    average_params(neuriteID, t, angle, parent);
    UG_COND_THROW(neuriteID >= m_vNeurites.size(), "Requested neurite ID which is not present.");

    const Neurite& neurite = m_vNeurites[neuriteID];
    const std::vector<BranchingRegion>& vBR = neurite.vBR;

    // vector for new position
    vector3 pos(0.0);

    // FOUR CASES can occur:
    // 1. We are at a branching point.
    // 2. We are well inside a regular piece of neurite.
    // 3. We are at the tip of a neurite.
    // 4. We are at the soma.

    // case 1: branching point?
    BranchingRegion cmpBR(t);
    std::vector<BranchingRegion>::const_iterator it =
        std::upper_bound(vBR.begin(), vBR.end(), cmpBR, CompareBranchingRegionEnds());

    if (it != vBR.end() && it->tstart < t)
        pos_on_surface_bp(pos, neurite, t, angle, it, parent, this);

    // case 2: normal neurite position
    else if (t < 1.0)
        pos_on_surface_neurite(pos, neurite, neuriteID, t, angle);

    // case 3: tip of neurite
    else
        pos_on_surface_tip(pos, neurite, parent, this);

    // case 4: soma
    // TODO: implement!


    // save new surface params for new vertex
    m_aaSurfParams[vrt].neuriteID = neuriteID;
    m_aaSurfParams[vrt].axial = t;
    m_aaSurfParams[vrt].angular = angle;

    // set position
    set_pos(vrt, pos);

    return 1.0;
}

// -------------------------------------------------------- //
// DO NOT CHANGE below this line! Needed for serialization. //

std::ostream& operator<<(std::ostream &os, const NeuriteProjector::SurfaceParams& surfParams)
{
    using std::ostringstream;
    ostringstream strs;
    strs << surfParams.neuriteID << " ";
    strs << surfParams.axial << " ";
    strs << surfParams.angular;
    os << strs.str();
    return os;
}

std::istream& operator>>(std::istream& in, NeuriteProjector::SurfaceParams& surfParams)
{
    std::string temp;
    using boost::lexical_cast;

    in >> temp;
    surfParams.neuriteID = lexical_cast<uint32>(temp);
    temp.clear();

    in >> temp;
    surfParams.axial = (lexical_cast<float>(temp));
    temp.clear();

    // onset
    in >> temp;
    surfParams.angular = lexical_cast<float>(temp);
    temp.clear();

    return in;
}


} // namespace ug
