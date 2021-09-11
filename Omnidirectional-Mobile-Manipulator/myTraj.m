function [T,Ve,p,pd,Q,rpy,rpyd,s,sd,sdd,c,cd,cdd] = myTraj(pstart,pfinish,Qstart,Qfinish,t,ts,qd0,qd1)
% Helper function for SLERP and cartesian trajectory planning
    [s,sd,sdd,c,cd,cdd] = tpolyMod(0,1,t,qd0,qd1);
    dp = pfinish-pstart;
    for i = 1:length(s)
        p(:,i) = pstart + s(i)*dp;
        if i < length(s)
            Q(:,i) = Qstart.interp(Qfinish,s(i));
        else
            Q(:,i) = Qstart.interp(Qfinish,1);
        end
%         pd(:,i) = (dp)*sd(i);
        pd(:,i) = (dp/norm(dp))*sd(i);
    end
    p = p';
    pd = pd';

    rpy = Q.torpy;
    rpyd(1,:) = zeros(1,3);
    for i = 1:length(s)-1
        if i < length(s)-1
            rpyd(i+1,:) = (rpy(i+1,:)-rpy(i,:))/ts;
        else
            rpyd(i+1,:) = zeros(1,3);
        end
    end
    

    T = rt2tr(rpy2r(rpy(:,:)),p(:,:)');
    Ve = [pd rpyd]';
end