def cleanSpecial(Object s) {
    // Remove bash special characters as defined in
    // https://www.howtogeek.com/439199/15-special-characters-you-need-to-know-for-bash/
    // and replace them with underscores
    return s.replace('~', '_')
        .replace('/', '_')
        .replace('#', '_')
        .replace('?', '_')
        .replace('*', '_')
        .replace('[]', '_')
        .replace('[', '_')
        .replace(']', '_')
        .replace(';', '_')
        .replace('&', '_')
        .replace('<', '_')
        .replace('>', '_')
        .replace('|', '_')
        .replace('!', '_')
        .replace('$', '_')
        .replace('{}', '_')
        .replace('{', '_')
        .replace('}', '_')
}
