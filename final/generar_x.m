function [ x, stateSeq ] = generar_x( hmm, C )
    ok = false;

    while ~ok
        [x,stateSeq] = genhmm(hmm);

        ok = true;

        for i=2:4
            if sum(stateSeq==i) < C
                ok = false;
                break;
            end

        end

    end
end

