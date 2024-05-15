%% Bistable map with chaos
%
% May 2024
% Peter Ashwin, Julian Newman, Raphael Romer
%
clear all;
close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% save pdfs
printfigs=true;

% which to run
runfigs0=true;
runfigs1=true;
runfigs2=true;
runfigs3=true;
runfigs4=true;
runfigs5=true;
runfigs6=true;
%runjnfigs=true;

% map parameter
% tangencies at beta=+-0.5463=tan(+-0.5) for alpha=2.
alpha=2;

% critical value of beta
betac=pi/2-1;

% fibre map
fx=@(x,b,s) alpha*atan(x)+b+s;

% fiber derivative
dfdx=@(x,b,s) alpha/(1.+x^2);

% inverse for given b and s
fxinv=@(x,b,s) tan((x-b-s)/alpha);

% noise map
fth=@(th) mod(3*th,2*pi);


if runfigs0==true
    %% plot map for different values of beta,noise
    f1=figure();
    fig=f1.Number;
    f1.Position(3:4)=[300 250];
    
    xx=-7:0.1:7;
    bs=[0.0 0.5 1.0];
    plot(xx,fx(xx,bs(1),0));
    hold on
    plot(xx,fx(xx,bs(2),0));
    plot(xx,fx(xx,bs(3),0))
    plot(xx,xx,'k');
    plot(xx,0,'k');
    xlabel("$x$");
    ylabel("$f(x)$");
    legend("$\beta$="+bs(1),"$\beta$="+bs(2),"$\beta$="+bs(3),"Location","se")
    hold off
    xlim([-4 4]);
    ylim([-4 4]);
    
    savefigure("figplot-v13",fig,printfigs)
end

if runfigs1==true
    %% plot example of attractors in phase space

    % noise amplitude
    sigma=0.24;
    %sigma=0.3;


    for beta=[0 0.1 0.2 0.3 0.33079 0.5]

        % plot time
        Nplot=5e6;

        % simulation time
        N=Nplot*100;

        % transient time - don't plot
        Ntrans=100;
        
        
        % initial condition
        x0=-2.33;
        th0=rand(1)*2*pi;
        
        ths=zeros(N,1);

        for i=1:N
            ths(i)=th0;
            th1=fth(th0);
            th0=th1;
        end

        xs=zeros(Nplot,2);
        xsc=zeros(N,1);

        %% find upper and lower attractors
        %
        for j=1:2
            x0s=[-2.3 2.3];
            x0=x0s(j);
            for i=1:Nplot
                xs(i,j)=x0;
                x1=fx(x0,beta,sigma*cos(ths(i)));
                x0=x1;
            end
        end

        %% find chaotic saddle over generic point
        x1=0.0;
        % compute backwards cocycle over computed *long* trajectory for theta
        for i=N-1:-1:1
            xsc(i+1)=x1;
            x0=fxinv(x1,beta,sigma*cos(ths(i)));
            x1=x0;
        end

        %% compute UPOs up to period Np
        % Np=9;
        % % ternary expansions up to period Np
        % upocode=table2array(combinations([0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2]));
        % Nupo=size(upocode,1);
        % thupo=zeros(Nupo,1);
        % for j=1:Nupo
        %     temp=0;
        %     for k=1:Np
        %         temp=temp+upocode(j,k)*3^(-k);
        %     end
        %     thupo(j)=2*pi*temp/(1-3^(-Np));
        % end
        % thupoinv=mod(3*thupo,2*pi);
        % 
        % xsthupo=zeros(Nupo,1);
        % % transient number of circuits to take
        % Nrep=10;
        % for kk=1:Nupo
        %     x1=0.0;
        %     % start with image at next point
        %     th1=thupo(kk);
        %     % no of steps backwards to end up back at same point as thupo(kk)
        %     for j=Nrep*Np+1:-1:1
        %         % one inverse step
        %         x0=fxinv(x1,beta,sigma*cos(th1));
        %         th0=fthupoinv(th1,thupo,thupoinv);
        %         x1=x0;
        %         th1=th0;
        %     end
        %    xsthupo(kk)=x1;
        % end

        %% plot attractors and saddle
        f1=figure();
        fig=f1.Number;
        f1.Position(3:4)=[400 350];

        % plot the chaotic saddle
        plot(xsc(Ntrans:Nplot-Ntrans),ths(Ntrans:Nplot-Ntrans),'.');
        hold on;
        % plot chaotic saddle for UPOs
        %plot(xsthupo,thupo,'.k');
        % plot one attractor
        plot(xs(Ntrans:Nplot,1),ths(Ntrans:Nplot),'.');
        % plot the other attractor
        plot(xs(Ntrans:Nplot,2),ths(Ntrans:Nplot),'.');
        
        xlim([-4 4]);
        ylim([0 2*pi]);
        xlabel("$x$");
        ylabel("$\theta$");
        title("$\beta=$"+beta+" $\sigma=$"+sigma);

        savefigure("figplotx-v13",fig,printfigs)

        % find LE on saddle
        le0=0;
        for i=1:N-1
            le0=le0+log(dfdx(xsc(i),beta,sigma*cos(ths(i))));
        end
        sprintf('fibre LE on saddle for beta=%d sigma=%d is %f',beta,sigma,(le0/N))

        if beta+sigma<=pi/2-1
            %% estimate box dimension of chaotic saddle if \beta+\sigma small enough
            % number of dyads
            P=10;
            % box for analysis centred on mean in small interval
            thmin=0.0;
            thmax=0.01;
            xav=mean(xsc(ths<0.02));
            xmin=xav-0.01;
            xmax=xav+0.01;
            
            ngx=2^P; ngy=2^P;
            grid=zeros(ngx,ngy);
            gridc=zeros(ngx,ngy);

            rx=linspace(xmin,xmax,ngx);
            ry=linspace(thmin,thmax,ngy);
            % boxes hit by typical orbit
            for i=Ntrans:N-Ntrans
                ix=floor((xsc(i)-rx(1))/(rx(end)-rx(1))*ngx)+1;
                iy=floor((ths(i)-ry(1))/(ry(end)-ry(1))*ngy)+1;
                if ix>0 && ix<=ngx && iy>0 && iy<=ngy
                    % number of hits
                    grid(ix,iy)=grid(ix,iy)+1;
                    % any hit
                    gridc(ix,iy)=1;
                end
            end
            % % boxes hit by UPO
            % for i=1:Nupo
            %     ix=floor((xsthupo(i)-rx(1))/(rx(end)-rx(1))*ngx)+1;
            %     iy=floor((thupo(i)-ry(1))/(ry(end)-ry(1))*ngy)+1;
            %     if ix>0 && ix<=ngx && iy>0 && iy<=ngy
            %         % any hit
            %         gridc(ix,iy)=1;
            %     end
            % end


            f1a=figure();
            fig=f1a.Number;
            f1a.Position(3:4)=[300 600];
            subplot(2,1,1)
            p1=pcolor(rx,ry,gridc');
            p1.EdgeColor="none";
            xlabel("$x$");
            ylabel("$\theta$");          
            title("$\beta=$"+beta+" $\sigma=$"+sigma);


            Pj=3;
            % Now count the boxes using sizes 2^-P..2^(Pj-P)
            epsj=2.^(0:1:Pj)*range(rx)/2^(P);
            Nj=zeros(size(epsj));
            Nj(1)=sum(sum(gridc));
            gridold=gridc;
            % count hit boxes in grid of size 2^(P-j)
            for j=1:Pj
                ng=2^(P-j);

                gridnew=zeros(ng,ng);
                for k1=1:ng
                    for k2=1:ng
                        % 2 by 2 grid to reduce
                        tempm=gridold(2*k1-1:2*k1,2*k2-1:2*k2);
                        gridnew(k1,k2)=max(max(tempm));
                    end
                end
                gridold=gridnew;
                Nj(j+1)=sum(sum(gridnew));
            end

            lm1=fitlm(log10(1./epsj),log10(Nj));
            slm1=lm1.Coefficients{:,:};
            sprintf("beta=%f sigma=%f d_b=%f SE %f",beta,sigma,slm1(2,1),slm1(2,2))

            subplot(2,1,2)
            plot(log10(1./epsj),log10(Nj),"-*");
            xlabel("$\log(1/\epsilon)$");
            ylabel("$\log(N(\epsilon))$");
            title("slope="+slm1(2,1));
            
            savefigure("figsaddle-dim-v13",fig,printfigs)


        end


    end


    %keyboard

    N=1000;
    xs=zeros(N,1);
    ts=zeros(N,1);
    bs=zeros(N,1);
    ths=zeros(N,1);

    %% plot figures with various non-zero sigma
    beta=0.3;
    sumf=figure();
    summaryfig=sumf.Number;
    sumf.Position(3:4)=[600 300];
    sigmarange=[0.0 0.25 0.5 0.75];
    epsilon=0.001;
    bn=@(n) epsilon*n;

    for sigma=sigmarange
        x0=-2.3;
        for i=1:N
            x1=fx(x0,bn(i),sigma*cos(th0));
            th1=fth(th0);
            xs(i)=x1;
            ts(i)=i;
            bs(i)=bn(i);
            ths(i)=th1;
            x0=x1;
            th0=th1;
        end

        figure(summaryfig);
        plot(bs,xs);
        hold on
        xlabel("$\beta$");
        ylabel("$x$");
    end
    xx=-4:0.05:4;
    betax=xx-alpha* atan(xx);
    plot(betax,xx);
    xlim([0 1]);

    hold off
    legend("$\sigma$="+sigmarange(1),"$\sigma$="+sigmarange(2),...
        "$\sigma$="+sigmarange(3),"$\sigma$="+sigmarange(4),"$\sigma=\epsilon=0$","Location","se");

    savefigure("figsummary-v13",summaryfig,printfigs);


end

%% plot parameter space vs size of attractor on fibre

if runfigs2==true

    parf=figure();
    parfig=gcf().Number;clf;
    parf.Position(3:4)=[700 220];

    betas=-1:0.02:1.;
    sigmas=0:0.02:2.;
    Nt=2000;
    N =50000;

    [bb,ss]=meshgrid(betas,sigmas);

    xpN=zeros(size(bb));
    xmN=zeros(size(bb));
    leN=zeros(size(bb));

    for i=1:length(betas)
        tic
        beta=betas(i);
        for j=1:length(sigmas)
            sigma=sigmas(j);
            theta0=rand*2*pi;
            xp0=3;
            xm0=-3;
            lep0=0;
            lem0=0;
            for k=1:N
                xp1=fx(xp0,beta,sigma*cos(th0));
                xm1=fx(xm0,beta,sigma*cos(th0));
                % compute Lyapunov exponent after transient
                if k>Nt
                    lep0=lep0+log(dfdx(xp0,beta,sigma*cos(th0)));
                    lem0=lem0+log(dfdx(xm0,beta,sigma*cos(th0)));
                end
                th1=fth(th0);
                th0=th1;
                xp0=xp1;
                xm0=xm1;
            end
            xpN(j,i)=xp0;
            xmN(j,i)=xm0;
            leN(j,i)=max(lep0/(N-Nt),lem0/(N-Nt));
        end
        toc
    end
    xd=xpN-xmN;
    subplot(1,2,1)
    % number of attractors found
    p1=pcolor(bb,ss,real(xd>1e-10)+1);
    p1.EdgeColor='none';
    xlabel('$\beta$')
    ylabel('$\sigma$');
    colorbar;
    xd=xpN-xmN;
    subplot(1,2,2)
    % most positive LE on the attractors
    p2=pcolor(bb,ss,leN);
    p2.EdgeColor='none';
    xlabel('$\beta$')
    ylabel('$\sigma$');
    colorbar;

    savefigure("figparams-v13",parfig,printfigs);
end



%% plot with UPOS for fixed sigma, beta and different ramp rates

if runfigs3==true

    epsilons=[0.0001 0.001 0.01];

    for epsilon=epsilons
        tic
        upof=figure();
        upofig=gcf().Number;clf;
        upof.Position(3:4)=[600 300];

        sigma=0.2;
        N=ceil(1/epsilon);
        %alpha=2;
        bn=@(n) epsilon*n;

        % UPO period to test
        Np=6;

        % ternary expansions up to period Np
        upocode=table2array(combinations([0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2]));
        Nupo=size(upocode,1);
        thupo=zeros(Nupo,1);
        for j=1:Nupo
            temp=0;
            for k=1:Np
                temp=temp+upocode(j,k)*3^(-k);
            end
            thupo(j)=2*pi*temp/(1-3^(-Np));
        end

        xs=zeros(N,1);
        ts=zeros(N,1);
        bs=zeros(N,1);
        ths=zeros(N,1);

        %%

        for j=1:Nupo
            th0=thupo(j);
            x0=-3.0;
            for i=1:N
                x1=fx(x0,bn(i),sigma*cos(th0));
                %% constrain to UPO
                th1=fthupo(th0,thupo);
                xs(i)=x1;
                ts(i)=i;
                bs(i)=bn(i);
                ths(i)=th1;
                x0=x1;
                th0=th1;
            end
            plot(bs,xs);
            hold on

        end

        % plot "typical" trajectory on top
        x0=-2.3;
        th0=rand*2*pi;
        for i=1:N
            x1=fx(x0,bn(i),sigma*cos(th0));
            th1=fth(th0);
            xs(i)=x1;
            ts(i)=i;
            bs(i)=bn(i);
            ths(i)=th1;
            x0=x1;
            th0=th1;
        end
        plot(bs,xs,'k',"LineWidth",2);

        % plot hysteresis curve on top
        xx=-4:0.05:4;
        betax=xx-alpha* atan(xx);
        plot(betax,xx,'k');
        plot(betax+sigma,xx,'k:');
        plot(betax-sigma,xx,'k:');

        xlabel("$\beta$");
        ylabel("$x$");
        xlim([0 1]);
        title("$\sigma=$"+sigma+", $\epsilon=$"+epsilon+", $\alpha=$"+alpha);
        hold off

        savefigure("figupo-v13",upofig,printfigs);
        toc
    end
end

%% plot escape times with UPOS for fixed sigma, beta and different ramp rates

if runfigs4==1
    % epsilons to be considered
    epsilons=[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2];

    Neps=length(epsilons);

    % ensemble size
    Nens=300;

    % UPO periods to test
    Np=2;

    % ternary expansions up to period Np
    %upocode=table2array(combinations([0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2]));
    upocode=table2array(combinations([0 1 2],[0 1 2]));
    Nupo=size(upocode,1);
    thupo=zeros(Nupo,1);
    for j=1:Nupo
        temp=0;
        for l=1:Np
            temp=temp+upocode(j,l)*3^(-l);
        end
        thupo(j)=2*pi*temp/(1-3^(-Np));
    end

    betadist=zeros(Neps,Nens);

    for k=1:Neps
        tic
        epsilon=epsilons(k);


        sigma=0.2;
        N=ceil(1/epsilon);
        alpha=2;
        bn=@(n) epsilon*n;

        xs=zeros(N,1);
        ts=zeros(N,1);
        bs=zeros(N,1);
        ths=zeros(N,1);

        %%

        for j=1:Nupo
            th0=thupo(j);
            x0=-2.3;
            xthresh=1.0;
            i=1;
            while x0<xthresh && i<N
                x1=fx(x0,bn(i),sigma*cos(th0));
                % constrain to UPO
                th1=fthupo(th0,thupo);
                xs(i)=x1;
                ts(i)=i;
                bs(i)=bn(i);
                ths(i)=th1;
                x0=x1;
                th0=th1;
                i=i+1;
            end
            imax=i-1;
        end

        % now plot trajectories on top
        for j=1:Nens
            x0=-2.3;
            xthresh=1;
            th0=rand*2*pi;
            i=1;
            while x0<xthresh && i<N
                x1=fx(x0,bn(i),sigma*cos(th0));
                th1=fth(th0);
                xs(i)=x1;
                ts(i)=i;
                bs(i)=bn(i);
                ths(i)=th1;
                x0=x1;
                th0=th1;
                i=i+1;
            end
            imax=i-1;
            betadist(k,j)=bs(imax);
        end
        toc
    end

    %% plot figure of cdfs
    cdff=figure();
    cdffig=gcf().Number;clf;
    cdff.Position(3:4)=[500 500];

    for j=1:Neps
        subplot(ceil(Neps/2),2,j)
        cdfplot(betadist(j,:));
        hold on
        plot([0.5708,0.5708]-sigma,[0,1],'-');
        plot([0.5708,0.5708],[0,1],'-');
        plot([0.5708,0.5708]+sigma,[0,1],'-');
        title("$\epsilon=$"+epsilons(j));
        xlim([0.3 0.8]);
        xlabel("$\beta$");
        ylabel("cdf");
    end

    %%

    savefigure("figcdfs-v13",cdffig,printfigs);
end

%% plot dynamic tipping windows for changing epsilon

if runfigs5==1
    % epsilons to be considered
    epsilons=[1e-6 1e-5 1e-4 5e-4 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3 8e-3 9e-3 1e-2];
    
    sigma=0.2;
        
    Neps=length(epsilons);

    % ensemble size
    Nens=300;

    % UPO periods to test
    Np=5;

    % ternary expansions up to period Np
    upocode=table2array(combinations([0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2]));
    %upocode=table2array(combinations([0 1 2],[0 1 2]));
    
    Nupo=size(upocode,1);
    thupo=zeros(Nupo,1);

    for j=1:Nupo
        temp=0;
        for l=1:Np
            temp=temp+upocode(j,l)*3^(-l);
        end
        thupo(j)=2*pi*temp/(1-3^(-Np));
    end

    % record tipping beta
    betaupo=zeros(Neps,Nupo);
    betadist=zeros(Neps,Nens);
    betaLquartile=zeros(Neps,1);
    betaUquartile=zeros(Neps,1);
    betamedian=zeros(Neps,1);
    
    for k=1:Neps
        tic
        epsilon=epsilons(k);

        N=ceil(1/epsilon);
        alpha=2;
        bn=@(n) epsilon*n;
        
        xs=zeros(N,1);
        ts=zeros(N,1);
        bs=zeros(N,1);
        ths=zeros(N,1);

        %% ramp with UPO forcing
        for j=1:Nupo
            th0=thupo(j);
            x0=-2.3;
            xthresh=1.0;
            i=1;
            while x0<xthresh && i<N
                x1=fx(x0,bn(i),sigma*cos(th0));
                % constrain to UPO
                th1=fthupo(th0,thupo);
                xs(i)=x1;
                ts(i)=i;
                bs(i)=bn(i);
                ths(i)=th1;
                x0=x1;
                th0=th1;
                i=i+1;
            end
            imax=i-1;
            betaupo(k,j)=bs(imax);
        end

        %% ramp with general forcing
        for j=1:Nens
            x0=-2.3;
            xthresh=1;
            th0=rand*2*pi;
            i=1;
            while x0<xthresh && i<N
                x1=fx(x0,bn(i),sigma*cos(th0));
                th1=fth(th0);
                xs(i)=x1;
                ts(i)=i;
                bs(i)=bn(i);
                ths(i)=th1;
                x0=x1;
                th0=th1;
                i=i+1;
            end
            imax=i-1;
            betadist(k,j)=bs(imax);
        end
        betamedian(k)=median(betadist(k,:));
        betaLquartile(k)=quantile(betadist(k,:),0.25);
        betaUquartile(k)=quantile(betadist(k,:),0.75);
        toc
    end


    %% plot figure dynamic tipping windows

    winf=figure();
    winfig=gcf().Number;clf;
    winf.Position(3:4)=[500 300];

    plot(betaupo,epsilons);
    hold on
    plot(betamedian,epsilons,'*k');
    plot(betamedian,epsilons,'k');
    plot(betaLquartile,epsilons,'*k');
    plot(betaLquartile,epsilons,'k');
    plot(betaUquartile,epsilons,'*k');
    plot(betaUquartile,epsilons,'k');

    plot(betac-sigma,0,'*');
    plot(betac,0,'*');
    plot(betac+sigma,0,'*');
    
    title("$\sigma=$"+sigma);
    %xlim([0.3 0.8]);
    xlabel("$\beta$");
    ylabel("$\epsilon$");

    %%

    savefigure("figwin-v13",winfig,printfigs);

    %keyboard

end

if runfigs6==1
    %% find escape time distributions with UPOS
    % epsilons & sigmas considered
    epsilons=[1e-6 1e-5 0.0001 0.001];
    sigmas=[0.006 0.019 0.06 0.195];
    %epsilons=[1e-6 1e-5 0.0001 0.01];
    %sigmas=[0.006 0.019 0.06 0.5];

    Neps=length(epsilons);

    % ensemble size
    Nens=300;

    betadist=zeros(Neps,Nens);

    for k=1:Neps
        tic
        epsilon=epsilons(k);
        sigma=sigmas(k);

        N=ceil(1/epsilon);
        alpha=2;
        bn=@(n) epsilon*n;

        % UPO period to test
        Np=2;

        % ternary expansions up to period Np
        %upocode=table2array(combinations([0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2]));
        upocode=table2array(combinations([0 1 2],[0 1 2]));
        Nupo=size(upocode,1);
        thupo=zeros(Nupo,1);
        for j=1:Nupo
            temp=0;
            for l=1:Np
                temp=temp+upocode(j,l)*3^(-l);
            end
            thupo(j)=2*pi*temp/(1-3^(-Np));
        end

        xs=zeros(N,1);
        ts=zeros(N,1);
        bs=zeros(N,1);
        ths=zeros(N,1);

        for j=1:Nupo
            th0=thupo(j);
            x0=-2.3;
            xthresh=1.0;
            i=1;
            while x0<xthresh && i<N
                x1=fx(x0,bn(i),sigma*cos(th0));
                % constrain to UPO
                th1=fthupo(th0,thupo);
                xs(i)=x1;
                ts(i)=i;
                bs(i)=bn(i);
                ths(i)=th1;
                x0=x1;
                th0=th1;
                i=i+1;
            end
            imax=i-1;
        end

        % now plot trajectories on top
        for j=1:Nens
            x0=-2.3;
            xthresh=1;
            th0=rand*2*pi;
            i=1;
            while x0<xthresh && i<N
                x1=fx(x0,bn(i),sigma*cos(th0));
                th1=fth(th0);
                xs(i)=x1;
                ts(i)=i;
                bs(i)=bn(i);
                ths(i)=th1;
                x0=x1;
                th0=th1;
                i=i+1;
            end
            imax=i-1;
            betadist(k,j)=bs(imax);
            %       plot(bs(1:imax),xs(1:imax),'k',"LineWidth",2);
        end

        % plot hysteresis curve on top of that
        xx=-4:0.05:4;
        betax=xx-alpha* atan(xx);
        toc
    end

    %% plot figure of cdfs
    cdff=figure();
    cdffig=gcf().Number;clf;
    cdff.Position(3:4)=[500 500];

    for j=1:Neps
        subplot(ceil(Neps/2),2,j)
        cdfplot(betadist(j,:));
        hold on
        betamin=betac-sigmas(j);
        betamax=betac+sigmas(j);
        plot([betamin,betamin],[0,1],'-');
        plot([betac,betac],[0,1],'-');
        plot([betamax,betamax],[0,1],'-');
        title("$\epsilon=$"+epsilons(j)+" $\sigma=$"+sigmas(j));
        xlim([betamin,betamax]);
        xlabel("$\beta$");
        ylabel("cdf");
    end

    savefigure("figcdfs2-v13",cdffig,printfigs);
    fitlm(log10(epsilons),log10(sigmas))

    llf=figure();
    llfig=gcf().Number;clf;
    llf.Position(3:4)=[400 350];
    plot(log10(sigmas),log10(epsilons),'*')
    hold on;
    xlabel("$\log(\sigma)$");
    ylabel("$\log(\epsilon)$");
    plot(log10(sigmas),1.986*log10(sigmas)-1.5831)
    legend("ensemble approx.","$\log(\epsilon)=1.986\log(\sigma)-1.5831$","Location","northwest")

    % epsilon is approx sigma^2 for median at betac
    savefigure("figll-v13",llfig,printfigs);

end

% if runjnfigs==true
%     %% find escape time distributions with UPOS
%     % epsilons & sigmas considered
%     epsilons=[2.5e-6 1e-5];
%     sigmas=[5e-3 1.4142e-2];
%     %epsilons=[1e-6 1e-5 0.0001 0.01];
%     %sigmas=[0.006 0.019 0.06 0.5];
% 
%     Neps=length(epsilons);
% 
%     % ensemble size
%     Nens=300;
% 
%     betadist=zeros(Neps,Nens);
% 
%     for k=1:Neps
%         tic
%         epsilon=epsilons(k);
%         sigma=sigmas(k);
% 
%         N=ceil(1/epsilon);
%         a=2;
%         bn=@(n) epsilon*n;
% 
%         % UPO period to test
%         Np=2;
% 
%         % ternary expansions up to period Np
%         %upocode=table2array(combinations([0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2],[0 1 2]));
%         upocode=table2array(combinations([0 1 2],[0 1 2]));
%         Nupo=size(upocode,1);
%         thupo=zeros(Nupo,1);
%         for j=1:Nupo
%             temp=0;
%             for l=1:Np
%                 temp=temp+upocode(j,l)*3^(-l);
%             end
%             thupo(j)=2*pi*temp/(1-3^(-Np));
%         end
% 
%         xs=zeros(N,1);
%         ts=zeros(N,1);
%         bs=zeros(N,1);
%         ths=zeros(N,1);
% 
%         for j=1:Nupo
%             th0=thupo(j);
%             x0=-2.3;
%             xthresh=1.0;
%             i=1;
%             while x0<xthresh && i<N
%                 x1=fx(x0,bn(i),sigma*cos(th0));
%                 % constrain to UPO
%                 th1=fthupo(th0,thupo);
%                 xs(i)=x1;
%                 ts(i)=i;
%                 bs(i)=bn(i);
%                 ths(i)=th1;
%                 x0=x1;
%                 th0=th1;
%                 i=i+1;
%             end
%             imax=i-1;
%         end
% 
%         % now plot trajectories on top
%         for j=1:Nens
%             x0=-2.3;
%             xthresh=1;
%             th0=rand*2*pi;
%             i=1;
%             while x0<xthresh && i<N
%                 x1=fx(x0,bn(i),sigma*cos(th0));
%                 th1=fth(th0);
%                 xs(i)=x1;
%                 ts(i)=i;
%                 bs(i)=bn(i);
%                 ths(i)=th1;
%                 x0=x1;
%                 th0=th1;
%                 i=i+1;
%             end
%             imax=i-1;
%             betadist(k,j)=bs(imax);
%             %       plot(bs(1:imax),xs(1:imax),'k',"LineWidth",2);
%         end
% 
%         % plot hysteresis curve on top of that
%         xx=-4:0.05:4;
%         betax=xx-alpha* atan(xx);
%         toc
%     end
% 
%     %% plot figure of cdfs
%     cdff=figure();
%     cdffig=gcf().Number;clf;
%     cdff.Position(3:4)=[200 600];
% 
%     for j=1:Neps
%         subplot(Neps,1,j)
%         cdfplot(betadist(j,:));
%         hold on
%         betamin=betac-sigmas(j);
%         betamax=betac+sigmas(j);
%         plot([betamin,betamin],[0,1],'-');
%         plot([betac,betac],[0,1],'-');
%         plot([betamax,betamax],[0,1],'-');
%         title("$\epsilon=$"+epsilons(j)+" $\sigma=$"+sigmas(j));
%         %xlim([betamin,betamax]);
%         xlim([betac-0.01 betac+0.01]);
%         if j==2
%             xlim([betac-0.02 betac+0.02]);
%         end
% 
% 
%         xlabel("$\beta$");
%         ylabel("cdf");
%     end
% 
%     savefigure("figcdfs2-v13",cdffig,printfigs);
%     fitlm(log10(epsilons),log10(sigmas))
% 
%     llf=figure();
%     llfig=gcf().Number;clf;
%     llf.Position(3:4)=[400 350];
%     plot(log10(sigmas),log10(epsilons),'*')
%     hold on;
%     xlabel("$\log(\sigma)$");
%     ylabel("$\log(\epsilon)$");
%     plot(log10(sigmas),1.986*log10(sigmas)-1.5831)
%     legend("ensemble approx.","$\log(\epsilon)=1.986\log(\sigma)-1.5831$","Location","northwest")
% 
%     % epsilon is approx sigma^2 for median at betac
%     savefigure("figll-v13",llfig,printfigs);
% 
% end

keyboard

%% find UPO in array thupo closest to image of th

% upo noise map
function thnext=fthupo(th,thupo)
temp=mod(3*th,2*pi);
[c, idx]=min(abs(temp-thupo));
thnext=thupo(idx);
end


% upo noise map inverse: find closest in thupo that is mapped to th
function thprev=fthupoinv(th,thupo,thupoinv)
[c, idx]=min(abs(th-thupoinv));
thprev=thupo(idx);
end
