# QAM-Audio

MATLAB development of QAM system to transmit data over an audio channel.  The receiver is currently unable to detect or correct for constellation rotation errors of -90, 180, or 90 degrees.

**QAM4 Screenshot.png** (below) shows the simulated result of a simple configuration so that the principle of operation of QAM can be understood from the plots.

![](https://raw.githubusercontent.com/xioTechnologies/QAM-Audio/master/Images/QAM4%20Screenshot.png)


**QAM64 Screenshot.png** (below) shows the result of QAM256 being used to send and receive 96 kb/s through and audio loopback created by connecting the sound card's line in and line out sockets.

![](https://raw.githubusercontent.com/xioTechnologies/QAM-Audio/master/Images/QAM256%20Screenshot.png)

**Block Diagram.png** (below) is a block diagram of the complete system.

![](https://raw.githubusercontent.com/xioTechnologies/QAM-Audio/master/Images/Block%20Diagram.png)
