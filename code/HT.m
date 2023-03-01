function [A,B]=HT(curveId, xy, b,t)

    switch curveId
        case 'cilindroZ'
            [A,B]=CilParam(xy);
        case 'cilindroY'
            [A,B]=CilParam(xy);
        case 'cilindroX'
            [A,B]=CilParam(xy);
        case 'piano'
            [A,B]=Piano(xy);
        case 'pianoX'
            [A,B]=PianoX(xy);
        case 'pianoY'
            [A,B]=PianoY(xy);
        case 'pianoTrasv'
            [A,B]=pianoTrasv(xy);
        case 'mConv'
            [A,B]=mConv(xy, b);
        case 'Cono'
            [A,B]=ConeAlpha(xy,t);
        case 'ConoX'
            [A,B]=ConoParamX(xy,t);
        case 'ConoY'
            [A,B]=ConoParamY(xy,t);
        case 'Sphere'
            [A,B]=SferaParam(xy);
        case 'TorusZ'
            [A,B]=TorusParamZ(xy);
        case 'TorusY'
            [A,B]=TorusParamY(xy);
        case 'TorusX'
            [A,B]=TorusParamX(xy);
        case 'circle'
            [A,B]=circle(xy);
        otherwise
            disp('The given family of curves is not included in the atlas')
    end

end