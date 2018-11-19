function val = allm(v)

val = 0;

for i = 1:length(v)
    if v(i) ~= -1
        val = 1;
    end
end

end