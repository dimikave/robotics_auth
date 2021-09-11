function [J,Ja,qd,Ve,X,hx,hy,hz,k] = mytrajLoop(q,qdot,hx_prev,hy_prev,hz_prev,k,dx,dy,TMB,lwr)

    hx(1) = hx_prev(k);
    hy(1) = hy_prev(k);
    hz(1) = hz_prev(k);
    k = 1;
    for i = 1:length(q)
        if k == length(q)+1
            k = length(q)
        end
        T0M = rt2tr(rotz(q(k,3)),[q(k,1);q(k,2);0]);
        T0B = T0M*TMB;
        Jp = [eye(2) [-dx*cos(q(k,3))-dy*sin(q(k,3));dx*sin(q(k,3))-dy*cos(q(k,3))];zeros(3);0 0 1];
        Ja = lwr.jacob0(q(k,4:9));
        J(:,:,k) = [Jp Ja];
        qd = qdot(k,:)';
        Ve(:,k) = J(:,:,k)*qd;



        g = lwr.fkine(q(k,4:9));
        X(:,:,k) = T0B*g.T;
        hx(k) = X(1,4,k);
        hy(k) = X(2,4,k);
        hz(k) = X(3,4,k);
        k = k+1;
    end
    k = k-1;
end