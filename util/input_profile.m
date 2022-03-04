function profile=input_profile(time,val)
    profile.signals.dimensions = 1;
    profile.time = time;
    profile.signals.values = val;
end