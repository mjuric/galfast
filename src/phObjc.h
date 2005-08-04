#if !defined(PHOBJC_H)
#define PHOBJC_H
/*
 * An enum for types of object
 */
typedef enum {
   OBJ_UNK = 0, 
   OBJ_CR, 
   OBJ_DEFECT, 
   OBJ_GALAXY, 
   OBJ_GHOST, 
   OBJ_KNOWNOBJ,
   OBJ_STAR, 
   OBJ_TRAIL,
   OBJ_SKY,
   OBJ_NTYPE
} OBJ_TYPE;

/* here are the values the "flags" field can take, in bitwise fashion */
#define OBJECT1_CANONICAL_CENTER 0x1	/* used canonical, not local, centre
				   pragma typedef { */
#define OBJECT1_BRIGHT 0x2	/* detected by Bright Objects */
#define OBJECT1_EDGE 0x4	/* object is too close to edge of frame */
#define OBJECT1_BLENDED 0x8	/* object is/was blended */
#define OBJECT1_CHILD  0x10	/* object is a child */
#define OBJECT1_PEAKCENTER 0x20	/* given centre is position of peak pixel */
#define OBJECT1_NODEBLEND 0x40	/* no deblending attempted */
#define OBJECT1_NOPROFILE 0x80	/* too small to estimate a profile */
#define OBJECT1_NOPETRO 0x100	/* no Petrosian radius */
#define OBJECT1_MANYPETRO 0x200	/* more than one Petrosian radius */

#define OBJECT1_NOPETRO_BIG 0x400 /* no Petrosian radius as object is too big*/
#define OBJECT1_DEBLEND_TOO_MANY_PEAKS 0x800 /* too many peaks to deblend */
#define OBJECT1_CR 0x1000	/* contains a CR pixel */
#define OBJECT1_MANYR50 0x2000	/* more than one 50% radius */
#define OBJECT1_MANYR90 0x4000	/* more than one 90% radius */
#define OBJECT1_BAD_RADIAL 0x8000 /* some low S/N radial points */
#define OBJECT1_INCOMPLETE_PROFILE 0x10000 /* r_P includes off-frame pixels */
#define OBJECT1_INTERP 0x20000	/* object contains interpolated pixels */
#define OBJECT1_SATUR 0x40000	/* object contains saturated pixels */
#define OBJECT1_NOTCHECKED 0x80000 /* object contains NOTCHECKED pixels */

#define OBJECT1_SUBTRACTED 0x100000 /* object had wings subtracted */
#define OBJECT1_NOSTOKES 0x200000 /* object has no measured stokes params */
#define OBJECT1_BADSKY 0x400000  /* sky level is so bad that object is -ve */
#define OBJECT1_PETROFAINT 0x800000	/* >= 1 Petrosian radius too faint */
#define OBJECT1_TOO_LARGE 0x1000000	/* object is too large */
#define OBJECT1_DEBLENDED_AS_PSF 0x2000000 /* deblender treated obj as PSF */
#define OBJECT1_DEBLEND_PRUNED 0x4000000 /* deblender pruned peak list */
#define OBJECT1_ELLIPFAINT 0x8000000	/* Centre's fainter than desired
					   elliptical isophote */
#define OBJECT1_BINNED1 0x10000000 /* object was found in 1x1 binned image */
#define OBJECT1_BINNED2 0x20000000 /* object was found in 2x2 binned image */

#define OBJECT1_BINNED4 0x40000000 /* object was found in 4x4 binned image */
#define OBJECT1_MOVED 0x80000000 /* Object may have moved; largest flag value
				     pragma } OBJECT1_FLAGS */

#define OBJECT1_DETECTED (OBJECT1_BINNED1|OBJECT1_BINNED2|OBJECT1_BINNED4)
/*
 * here are values for flags2
 */
#define OBJECT2_DEBLENDED_AS_MOVING 0x1	/* deblended as a moving object
					   pragma typedef { */
#define OBJECT2_NODEBLEND_MOVING 0x2	/* no deblend of moving object */
#define OBJECT2_TOO_FEW_DETECTIONS 0x4  /* too few detections to deblend
					   as moving */
#define OBJECT2_BAD_MOVING_FIT 0x8	/* fit to moving object was too poor */
#define OBJECT2_STATIONARY 0x10		/* velocity is consistent with zero */
#define OBJECT2_PEAKS_TOO_CLOSE 0x20	/* at least some peaks were too
					   close, and thus merged */
#define OBJECT2_BINNED_CENTER 0x40	/* image was binned while centroiding*/
#define OBJECT2_LOCAL_EDGE 0x80		/* per-band centre's too near edge */
#define OBJECT2_BAD_COUNTS_ERROR 0x100	/* psf|fiberCountsErr is bad/unknown */
#define OBJECT2_BAD_MOVING_FIT_CHILD 0x200 /* moving child's fit was too poor*/


#define OBJECT2_DEBLEND_UNASSIGNED_FLUX 0x400 /* deblender failed to assign
						 enough of flux to children */
#define OBJECT2_SATUR_CENTER 0x800	/* object's centre's saturated */
#define OBJECT2_INTERP_CENTER 0x1000	/* object's centre's interpolated */
#define OBJECT2_DEBLENDED_AT_EDGE 0x2000 /* object's deblended although EDGE */
#define OBJECT2_DEBLEND_NOPEAK 0x4000	/* object had no detected peak */
#define OBJECT2_PSF_FLUX_INTERP 0x8000	/* a significant amount of
					   PSF's flux is interpolated */
#define OBJECT2_TOO_FEW_GOOD_DETECTIONS 0x10000	/* too few good detections to
						   deblend as moving */
#define OBJECT2_CENTER_OFF_AIMAGE 0x20000 /* at least one peak's centre lay off
					     the atlas image in some band */
#define OBJECT2_DEBLEND_DEGENERATE 0x40000 /* at least one potential child has
					      been pruned as being too similar
					      to some other template */


#define OBJECT2_BRIGHTEST_GALAXY_CHILD 0x80000 /* this is the brightest child
						  galaxy in a blend */
#define OBJECT2_CANONICAL_BAND 0x100000 /* This band was primary (usually r')*/
#define OBJECT2_AMOMENT_UNWEIGHTED 0x200000 /* `adaptive' moments are actually
					       unweighted */
#define OBJECT2_AMOMENT_SHIFT 0x400000	/* centre moved too far while
					   determining adaptive moments */
#define OBJECT2_AMOMENT_MAXITER 0x800000 /* Too many iterations while
					    determining adaptive moments */
#define OBJECT2_MAYBE_CR 0x1000000	/* object may be a cosmic ray */
#define OBJECT2_MAYBE_EGHOST 0x2000000	/* object may be an electronics ghost*/
#define OBJECT2_NOTCHECKED_CENTER 0x4000000 /* object's centre is NOTCHECKED */
#define OBJECT2_HAS_SATUR_DN 0x8000000	/* Counts include DN in bleed trails */
#define OBJECT2_SPARE4 0x10000000	/* unused */
#define OBJECT2_SPARE3 0x20000000	/* unused */
#define OBJECT2_SPARE2 0x40000000	/* unused */
#define OBJECT2_SPARE1 0x80000000	/* unused 
					   pragma } OBJECT2_FLAGS */
/*
 * These are book-keeping bits, and aren't written to disk
 */
#define OBJECT3_HAS_SATUR_DN 0x1	/* object has extra (saturated) DN
					   pragma typedef { */
#define OBJECT3_MEASURED 0x10000000	/* object has been measured */
#define OBJECT3_GROWN_MERGED 0x20000000	/* growing led to a merger */
#define OBJECT3_HAS_CENTER 0x40000000	/* OBJC has a canonical centre */
#define OBJECT3_MEASURE_BRIGHT 0x80000000 /* object should be measured bright
					     pragma } OBJECT3_FLAGS */

#endif

