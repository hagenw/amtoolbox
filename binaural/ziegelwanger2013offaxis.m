function y=ziegelwanger2013offaxis(p,x)
%ZIEGELWANGER2013OFFAXIS
%   y=ziegelwanger2013offaxis(p,x)
%
%   Input:
%       p: 
%           r.......sphere radius
%           xM......sphere offset on x-axis
%           yM......sphere offset on y-axis
%           zM......sphere offset on z-axis
%           delay...constant delay
%           phi.....position of the ear (azimuth angle)
%           theta...position of the ear (elevation angle)
%
%       x:
%
%   Output:
%       y:

% AUTHOR: Harald Ziegelwanger, Acoustics Research Institute, Vienna, Austria


    global ch

    if length(p)==2
        global r
        global yM
        global delay
        xM=p(1);
        zM=p(2);
    else if length(p)==3
            r=p(1);
            xM=0;
            yM=p(2);
            zM=0;
            delay=p(3);
        else
            r=p(1);
            xM=p(2);
            yM=p(3);
            zM=p(4);
            delay=p(5);
        end
    end
    
    if length(p)<6
        yo=sqrt(r^2-xM^2-zM^2)+yM;
        thetao=acos(sqrt(r^2-zM^2)/r);
    else
        thetao=p(7);
    end
    
    M=sqrt(xM^2+yM^2+zM^2);

if ch==1
    %linkes Ohr----------------------------------------
    if length(p)<6
        phil=abs(atan((yo-yM)/-xM));
        if xM>0
            phil=pi-phil;
        end
    else
        phil=p(6);
    end
    
    beta=acos(-cos(x(:,2)).*(xM*cos(x(:,1))+yM*sin(x(:,1)))-zM*sin(x(:,2)));
    s2=-r+M*cos(beta)+sqrt(r^2+M^2*cos(beta).^2+2*M*r);
    gamma=pi-beta-acos((2*M^2+2*M*r-2*r*s2-s2.^2)/(2*M^2+2*M*r));
    if M==0
        s1=zeros(size(x,1),1);
    else
        s1=M*cos(beta)./(2*(M+r).*tan(gamma/2));
    end
    
    y=1/343*((r* ...
        ((sign(sin(thetao).*sin(x(:,2))+cos(thetao).*cos(x(:,2)).*cos(phil-x(:,1)))/2+0.5).* ...
        (1-sin(thetao).*sin(x(:,2))-cos(thetao).*cos(x(:,2)).*cos(phil-x(:,1)))+ ...
        (-sign(sin(thetao).*sin(x(:,2))+cos(thetao).*cos(x(:,2)).*cos(phil-x(:,1)))/2+0.5).* ...
        (1+acos(sin(thetao).*sin(x(:,2))+cos(thetao)*cos(x(:,2)).*cos(phil-x(:,1)))-pi/2))) ...
        +s1+s2) ...
        +delay-(M+r)/343;
end

if ch==2
    %rechtes Ohr----------------------------------------
    if length(p)<6
        if xM<=0
            phir=2*pi-abs(atan((yo-yM)/-xM));
        else
            phir=pi+abs(atan((yo-yM)/-xM));
        end
    else
        phir=p(6);
    end
    
    beta=acos(-cos(x(:,2)).*(xM*cos(x(:,1))+yM*sin(x(:,1)))-zM*sin(x(:,2)));
    s2=-r+M*cos(beta)+sqrt(r^2+M^2*cos(beta).^2+2*M*r);
    gamma=pi-beta-acos((2*M^2+2*M*r-2*r*s2-s2.^2)/(2*M^2+2*M*r));
    if M==0
        s1=zeros(size(x,1),1);
    else
        s1=M*cos(beta)./(2*(M+r).*tan(gamma/2));
    end
    
    y=1/343*((r* ...
        ((sign(sin(thetao).*sin(x(:,2))+cos(thetao).*cos(x(:,2)).*cos(phir-x(:,1)))/2+0.5).* ...
        (1-sin(thetao).*sin(x(:,2))-cos(thetao).*cos(x(:,2)).*cos(phir-x(:,1)))+ ...
        (-sign(sin(thetao).*sin(x(:,2))+cos(thetao).*cos(x(:,2)).*cos(phir-x(:,1)))/2+0.5).* ...
        (1+acos(sin(thetao).*sin(x(:,2))+cos(thetao)*cos(x(:,2)).*cos(phir-x(:,1)))-pi/2))) ...
        +s1+s2) ...
        +delay-(M+r)/343;
end
end