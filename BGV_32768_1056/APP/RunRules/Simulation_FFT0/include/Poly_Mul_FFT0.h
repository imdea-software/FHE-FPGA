/**\file */
#ifndef SLIC_DECLARATIONS_Poly_Mul_FFT0_H
#define SLIC_DECLARATIONS_Poly_Mul_FFT0_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define Poly_Mul_FFT0_N (32768)
#define Poly_Mul_FFT0_TYPE_SIZEBYTES (16)
#define Poly_Mul_FFT0_n (15)
#define Poly_Mul_FFT0_PCIE_ALIGNMENT (16)


/*----------------------------------------------------------------------------*/
/*--------------------------- Interface writeLMem ----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'writeLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size param_size bytes.
 */
void Poly_Mul_FFT0_writeLMem(
	int32_t param_size,
	int32_t param_start,
	const int64_t *instream_fromcpu);

/**
 * \brief Basic static non-blocking function for the interface 'writeLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size param_size bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_writeLMem_nonblock(
	int32_t param_size,
	int32_t param_start,
	const int64_t *instream_fromcpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'writeLMem'
 * 
 */
typedef struct { 
	int32_t param_size; /**<  [in] Interface Parameter "size". */
	int32_t param_start; /**<  [in] Interface Parameter "start". */
	const int64_t *instream_fromcpu; /**<  [in] The stream should be of size param_size bytes. */
} Poly_Mul_FFT0_writeLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Poly_Mul_FFT0_writeLMem_run(
	max_engine_t *engine,
	Poly_Mul_FFT0_writeLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'writeLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_writeLMem_run_nonblock(
	max_engine_t *engine,
	Poly_Mul_FFT0_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Poly_Mul_FFT0_writeLMem_run_group(max_group_t *group, Poly_Mul_FFT0_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_writeLMem_run_group_nonblock(max_group_t *group, Poly_Mul_FFT0_writeLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Poly_Mul_FFT0_writeLMem_run_array(max_engarray_t *engarray, Poly_Mul_FFT0_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_writeLMem_run_array_nonblock(max_engarray_t *engarray, Poly_Mul_FFT0_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Poly_Mul_FFT0_writeLMem_convert(max_file_t *maxfile, Poly_Mul_FFT0_writeLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface readLMem ----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'readLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size param_size bytes.
 */
void Poly_Mul_FFT0_readLMem(
	int32_t param_size,
	int32_t param_start,
	int64_t *outstream_tocpu);

/**
 * \brief Basic static non-blocking function for the interface 'readLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size param_size bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_readLMem_nonblock(
	int32_t param_size,
	int32_t param_start,
	int64_t *outstream_tocpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'readLMem'
 * 
 */
typedef struct { 
	int32_t param_size; /**<  [in] Interface Parameter "size". */
	int32_t param_start; /**<  [in] Interface Parameter "start". */
	int64_t *outstream_tocpu; /**<  [out] The stream should be of size param_size bytes. */
} Poly_Mul_FFT0_readLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Poly_Mul_FFT0_readLMem_run(
	max_engine_t *engine,
	Poly_Mul_FFT0_readLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'readLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_readLMem_run_nonblock(
	max_engine_t *engine,
	Poly_Mul_FFT0_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Poly_Mul_FFT0_readLMem_run_group(max_group_t *group, Poly_Mul_FFT0_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_readLMem_run_group_nonblock(max_group_t *group, Poly_Mul_FFT0_readLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Poly_Mul_FFT0_readLMem_run_array(max_engarray_t *engarray, Poly_Mul_FFT0_readLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_readLMem_run_array_nonblock(max_engarray_t *engarray, Poly_Mul_FFT0_readLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Poly_Mul_FFT0_readLMem_convert(max_file_t *maxfile, Poly_Mul_FFT0_readLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] instream_inVector1 The stream should be of size 524288 bytes.
 * \param [in] instream_inVector2 The stream should be of size 524288 bytes.
 * \param [out] outstream_resultBNH The stream should be of size (524288 / 2) bytes.
 * \param [out] outstream_resultNHN The stream should be of size (524288 / 2) bytes.
 */
void Poly_Mul_FFT0(
	const uint64_t *instream_inVector1,
	const uint64_t *instream_inVector2,
	uint64_t *outstream_resultBNH,
	uint64_t *outstream_resultNHN);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] instream_inVector1 The stream should be of size 524288 bytes.
 * \param [in] instream_inVector2 The stream should be of size 524288 bytes.
 * \param [out] outstream_resultBNH The stream should be of size (524288 / 2) bytes.
 * \param [out] outstream_resultNHN The stream should be of size (524288 / 2) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_nonblock(
	const uint64_t *instream_inVector1,
	const uint64_t *instream_inVector2,
	uint64_t *outstream_resultBNH,
	uint64_t *outstream_resultNHN);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	const uint64_t *instream_inVector1; /**<  [in] The stream should be of size 524288 bytes. */
	const uint64_t *instream_inVector2; /**<  [in] The stream should be of size 524288 bytes. */
	uint64_t *outstream_resultBNH; /**<  [out] The stream should be of size (524288 / 2) bytes. */
	uint64_t *outstream_resultNHN; /**<  [out] The stream should be of size (524288 / 2) bytes. */
} Poly_Mul_FFT0_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Poly_Mul_FFT0_run(
	max_engine_t *engine,
	Poly_Mul_FFT0_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_run_nonblock(
	max_engine_t *engine,
	Poly_Mul_FFT0_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Poly_Mul_FFT0_run_group(max_group_t *group, Poly_Mul_FFT0_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_run_group_nonblock(max_group_t *group, Poly_Mul_FFT0_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Poly_Mul_FFT0_run_array(max_engarray_t *engarray, Poly_Mul_FFT0_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Poly_Mul_FFT0_run_array_nonblock(max_engarray_t *engarray, Poly_Mul_FFT0_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Poly_Mul_FFT0_convert(max_file_t *maxfile, Poly_Mul_FFT0_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* Poly_Mul_FFT0_init(void);

/* Error handling functions */
int Poly_Mul_FFT0_has_errors(void);
const char* Poly_Mul_FFT0_get_errors(void);
void Poly_Mul_FFT0_clear_errors(void);
/* Free statically allocated maxfile data */
void Poly_Mul_FFT0_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int Poly_Mul_FFT0_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int Poly_Mul_FFT0_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_Poly_Mul_FFT0_H */

