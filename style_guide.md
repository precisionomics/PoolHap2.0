# Java Style Guide:
To encourage consistent formatting, below details some conventions I hope to to implement moving forward.
Much of it is drawn from the [Twitter Java Style Guide](https://github.com/twitter/commons/blob/master/src/java/com/twitter/common/styleguide.md).

### Formatting
#### 100 column limit
A balance between fewer linebreaks and ease of reading.

#### No tab characters: use 4 spaces instead
Tab characters can cause a lot of issues. Instead, for indentations, 4 space characters should be used.
Use an IDE function to convert existing tab characters to 4 space characters and use 4 space characters for future indentations.

#### No trailing whitespace
Use an IDE function to automatically trim trailing whitespaces.

#### Indentation style
For consistency, use the "one true brace style" with 4 column indent sizes.
```java
// Correct.
if (true) {
    foo(x);
} else {
    bar(x);
}

// Don't do this.
if (true)
    foo(x);

// Or this.
if (false) bar(x);
```

Continuation indent is again at 4 columns. Nested continuations also at 4 columns.
```java
// Example from Twitter's styleguide.md

// Bad.
//   - Line breaks are arbitrary.
//   - Scanning the code makes it difficult to piece the message together.
throw new IllegalStateException("Failed to process request" + request.getId()
    + " for user " + user.getId() + " query: '" + query.getText()
    + "'");

// Good.
//   - Each component of the message is separate and self-contained.
//   - Adding or removing a component of the message requires minimal reformatting.
throw new IllegalStateException("Failed to process"
    + " request " + request.getId()
    + " for user " + user.getId()
    + " query: '" + query.getText() + "'");
    
```

For long method declarations:
```java
// Do this.
public String foo(
    int a,
    int b,
    boolean c,
    double d) {

    ...
}
```

For long chained method calls, wrap with leading periods at the 2nd dot and a linebreak at the elind:
```java
String a = ObjectA.methodB()
    .helperC()
    .propertyD()
    .somethingE();
    
```

#### Operators
Please pad them with spaces.
```java
// For example...
int a = b + c + 5;

// Don't do this
int a=b+c+5;
```

Be explicit about operation ordering, even if it's obvious.
```java
// Don't do this
int c = n + q * 5;

// Do this
int c = n + (q * 5);
```

When wrapping, same as chained method calls, wrap with leading operators.
```java
System.out.println("testestestestes"
    + a
    + b
    + c);
    
```

#### Comments
* Javadoc above methods with `@param`, `@throws`, and `@return` tags as necessary.
* Javadoc above methods with `@author` and `@since` tags as necessary.

In-line comments: 1 space character, no capital and no period.
```java
int a = b + c; // in-line
```

If in-line comments exceeds column limit, move to above the line with the same indentation
and proper punctuation.
```java
// For example, this comment.
int a = b + c;
```

If comments are across many lines and warrant a comment block, format as follows:
```java
// Don't do this
/* Something
something
something
something */

// Do this
/*
 *  2 white spaces after asterisk
 *  something
 *  something
 */
```

Overall, avoid unnecessary newlines. Use common sense for logical code chunking. **Prioritize
readability over formal conventions.**

### Variable Naming
PascalCase for classes/types, camelCase for variables, UPPER_SNAKE for constants.

**Important Note**: please avoid snake_case for Java. My brain keeps thinking I'm reading Python :(
    
#### Make clear and concise variable names
Avoid using abbreviations that do not make sense at first glance. Prioritize clarity over conscision.
If logical abbreviations are used, comments are still appreciated to tell us what it is.
```java
// Bad
String lrrhn = "Annie"

// Acceptable
double hapFreq = 0.5 // haplotype frequency

// Best
String bigBadWolfName = "Jerry?"
```
