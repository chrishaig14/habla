function seq = obtener_sec_modelos( stateSeq )
    seq = [];
    for i = 1:length(stateSeq)-1
        if stateSeq(i) == 4 && stateSeq(i+1) ~= 4
            seq(length(seq)+1) = 4;
        end
        if stateSeq(i) == 7 && stateSeq(i+1) ~= 7
            seq(length(seq)+1) = 6;
        end
    end
end

