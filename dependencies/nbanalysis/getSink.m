function [presink, postsink,sink,chans] = getSink(sessdate)

presink = [];
postsink = [];
sink = [];
switch sessdate
    case '170413'
        presink       = 16;
        postsink      = 6;
    case '170411'
        presink       = 15;
        postsink      = 7;
    case '170324'
        presink       = 21;
        postsink      = 16;
    case '160307'
        presink       = 15;
        postsink      = 14;
    case '160308'
        presink       = 19;
        postsink      = 19;
    case '160226'
        presink       = 19;
        postsink     = 19;
    case '140807'
        sink          = 17;
        chans         = [4:21];
    case '140724'
        sink          = 21;
        chans         = [7:24];
    case '140725'
        sink          = 20;
        chans         = [7:24];
    case '140818'
        sink          = 11;
        chans         = [2:17];
        
    case '140820'
        sink          = 15;
        chans         = [5:22];
    case '140731'
        sink          = 13;
        chans = [4:21];
    case '140729'
        sink          = 16;
        chans = [1:24];
        
end