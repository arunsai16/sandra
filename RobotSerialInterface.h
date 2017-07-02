/** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/
/**                                                             **/
/**  DO NOT MODIFY THIS FILE!  DOING SO MAY IMPACT YOUR         **/
/**    PROGRAM'S ABILITY TO COMMUNICATE WITH THE RV-M1 ROBOT    **/
/**    CONTROLLER!!                                             **/
/**                                                             **/
/** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! **/


#define	RVM1_WriteToFile	0x10
#define RVM1_WriteToScreen	0x20
#define RVM1_WriteToComm	0x80

/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_InitializeComm( unsigned int, char * ) is used to     **/
/**    initialize the communications port to the RV-M1 Robot    **/
/**    and the log file.  The user specifices whether to        **/
/**    communicate with the Robot, write to the log file, write **/
/**    to the screen, or any combination of the three options.  **/
/**    If writing to a log file is selected, the user can       **/
/**    provide the name of the log file.  If a path is          **/
/**    the directory must already exist otherwise the file will **/
/**    not be created.                                          **/
/**                                                             **/
/**    To select the output option, choose one or more of the   **/
/**    following (bit-wise OR for 2 or more options):           **/
/**         RVM1_WriteToFile    to write to a log file          **/
/**         RVM1_WriteToScreen	to write to the screen          **/
/**         RVM1_WriteToComm	to send commands to the RV-M1   **/
/**                             Robot Controller                **/
/**                                                             **/
/**    Example:                                                 **/
/**      unsigned int OutputOptions;                            **/
/**      unsigned int RobotResult;                              **/
/**      OutputOptions = RVM1_WriteToFile | RVM1_WriteToComm;   **/
/**      RobotResult = RVM1_InitializeComm( OutputOptions ,     **/
/**                          "c:\\RVM1_Logs\\RobotLog.txt" );   **/
/**                                                             **/
/**      Send command to both the RV-M1 Controller and to the   **/
/**        log file "RobotLog.txt" located in the directory     **/
/**        "C:\RVM1_Logs".  The directory must exist prior to   **/
/**        program execution and if a file with that name       **/
/**        already exists, it will be over-written.             **/
unsigned int RVM1_InitializeComm( unsigned int, char * );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_CleanUp( void ) is used at the end of the user's      **/
/**    robot program to clean-up any open files and the         **/
/**    communications port to the RV-M1 Robot Controller.       **/
void RVM1_CleanUp( void );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_FlushBuffers( void ) simply forces all open file and  **/
/**    communications port buffers to be flushed.               **/
void RVM1_FlushBuffers( void );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_GripOpen( void ) causes the gripper to open.          **/
void RVM1_GripOpen( void );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**  RVM1_GripClose( void ) causes the gripper to close.        **/
void RVM1_GripClose( void );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_MoveJoint( float, float, float, float, float ) is     **/
/**    used to command the joints of the robot to new rotations **/
/**    relative to the current position.  The joints, in order: **/
/**                                                             **/
/**         Waist                                               **/
/**         Shoulder                                            **/
/**         Elbow                                               **/
/**         Wrist Pitch                                         **/
/**         Wrist Rotation                                      **/
void RVM1_MoveJoint( float, float, float, float, float );


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_Nest( void ) commands the robot to the "nest"         **/
/**    position for storage.                                    **/
void RVM1_Nest( void );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_Origin( void ) commands the robot to the "origin"     **/
/**    position.  That is Waist = 0., Shoulder = 0.,            **/
/**    Elbow = 0., Wrist Pitch = 0., and Wrist Roatation = 0.   **/
void RVM1_Origin( void );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_Reset( void ) reset the RV-M1 Robot controller.       **/
void RVM1_Reset( void );

/** ----------------------------------------------------------- **/


/** ----------------------------------------------------------- **/
/**                                                             **/
/**  RVM1_Speed( int ) sets the movement speed of the robot.    **/
/**    The input is and integer in the range of 1 (slow) to     **/
/**    9 (fast).                                                **/
void RVM1_Speed( int );

/** ----------------------------------------------------------- **/






