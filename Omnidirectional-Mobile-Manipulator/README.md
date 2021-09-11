----------------------------------
# Robotics Course Semester Project
----------------------------------

## Description:

This project is about the trajectory planning of
an omni-directional mobile manipulator robot in a
pick-and-place task of an object on a table. The trajectories
are planned both in the Joint Space and in the 
Task Space with quintic polynomial paths.

## Code:
Set-up functions for the simulation:
- lwr_create.p		 : Creates a robot arm object from the Matlab's
			   Robotics Toolbox (RTB) of Peter Corke.
- DrawCuboid.m 		 : Draws a cuboid in the currend figure
			   using 8 rectangular faces to demonstrate
			   the platform and the table.
- plotCylinderWithCaps.m : Plots a cylinder with specified parameters
			   to demonstrate the target-object.
- Robot.plot		 : To demonstrate the robot arm given by lwr object

Helper functions for the trajectory planning and other calculations:
- jtrajMod.m   		 : Modified version of the (RTB) function jtraj
		 	   to return the coefficients of the planned trajectories
- tpolyMod.m   		 : Modified version of the tpoly function of RTB to return the
		 	   coefficients of the planned trajectory
- myTraj       		 : Helper function for interpolating Poses with SLERP
			   Used in MainTS.m
- mytrajLoop		 : Helper function for calculations and forward kinematics
			   Used in MainJS.m


Main files:
- MainJS.m     		 : Main file - Trajectory planning in the JOINT SPACE
		 	   for the project's pick and place task. (Not all requirements
			   fulfilled).

- Main.m		 : Main file - Trajectory planning in the TASK SPACE with
			   project's requirements fulfilled.
## Report:
9351_KavelidisFrantzis_Report.pdf : Report of the project
